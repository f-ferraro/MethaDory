#' PyTorch NNET inference for MethaDory.
#'
#' Mirrors pipelines/05c_vc_optimization_nnet_test.py: reconstruct each
#' FlexibleNNet checkpoint from its `model_config` and run sigmoid forward on
#' the matching probes. Uses reticulate to call torch from R.

# Lazy reticulate handle. Set up the Python side (torch + the FlexibleNNet
# class) only once per session.
.nnet_py <- new.env(parent = emptyenv())

#' Initialize Python (torch + FlexibleNNet definition). Idempotent.
init_nnet_python <- function() {
  if (!is.null(.nnet_py$torch)) return(invisible(TRUE))

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for NNET inference but is not installed.")
  }

  # Force reticulate to use the pixi environment's Python, robustly across
  # machines. Two failure modes to avoid:
  #   1. If RETICULATE_PYTHON is unset, recent reticulate auto-downloads a
  #      uv-managed interpreter; on Apple Silicon that's an arm64 build that is
  #      ABI-incompatible with the x86_64 pixi/R env (dlopen: "incompatible
  #      architecture").
  #   2. Pointing it straight at $CONDA_PREFIX/bin/python makes reticulate
  #      classify the env as conda (it has a conda-meta/ dir) and demand a
  #      `conda` binary to activate it -- which fails on machines without
  #      Anaconda/conda installed ("Unable to find conda binary").
  # The pixi env is already fully activated by `pixi run`, so we point
  # reticulate at a symlink to the interpreter in a plain directory (no
  # conda-meta sibling): same Python, but treated as a regular interpreter, so
  # no conda is ever needed.
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_prefix) && !reticulate::py_available(initialize = FALSE)) {
    real_py <- file.path(conda_prefix, "bin", "python")
    if (file.exists(real_py)) {
      target_py <- real_py
      shim_dir <- file.path(tempdir(), "methadory_pyshim", "bin")
      dir.create(shim_dir, recursive = TRUE, showWarnings = FALSE)
      shim_py <- file.path(shim_dir, "python")
      if (!file.exists(shim_py)) {
        if (isTRUE(tryCatch(file.symlink(real_py, shim_py),
                            error = function(e) FALSE))) {
          target_py <- shim_py
        }
      } else {
        target_py <- shim_py
      }
      Sys.setenv(RETICULATE_PYTHON = target_py)
    }
  }

  torch <- tryCatch(reticulate::import("torch", convert = FALSE),
                    error = function(e) stop(
                      "Could not import Python 'torch'. Install pytorch in the ",
                      "active reticulate environment. Original error: ", e$message))
  reticulate::import("torch.nn", convert = FALSE)
  reticulate::import("numpy", convert = FALSE)

  reticulate::py_run_string("
import torch
import torch.nn as nn
import numpy as np

class FlexibleNNet(nn.Module):
    def __init__(self, n_features, hidden_sizes, weight_decay=0.0):
        super().__init__()
        self.activation = nn.Sigmoid()
        layers = []
        prev = n_features
        for h in hidden_sizes:
            layers.append(nn.Linear(prev, h))
            prev = h
        self.hidden_layers = nn.ModuleList(layers)
        self.output = nn.Linear(prev, 1)

    def forward(self, x):
        h = x
        for layer in self.hidden_layers:
            h = self.activation(layer(h))
        return self.output(h).squeeze(-1)

    def predict_proba(self, x):
        with torch.no_grad():
            return torch.sigmoid(self.forward(x))

_DEVICE = torch.device('cpu')

def load_nnet_checkpoint(path):
    try:
        ckpt = torch.load(path, map_location=_DEVICE, weights_only=False)
    except TypeError:
        ckpt = torch.load(path, map_location=_DEVICE)
    cfg = ckpt['model_config']
    model = FlexibleNNet(
        n_features=cfg['n_features'],
        hidden_sizes=cfg['hidden_sizes'],
        weight_decay=cfg.get('weight_decay', 0.0),
    )
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    meta = ckpt.get('metadata', {}) or {}
    return {
        'model': model,
        'probe_ids': list(cfg.get('probe_ids') or []),
        'signature': meta.get('signature'),
        'fold_id': meta.get('fold_id', ''),
    }

def nnet_predict(model, X_np):
    X = torch.from_numpy(X_np.astype('float32'))
    return model.predict_proba(X).numpy()
")

  .nnet_py$torch <- torch
  .nnet_py$py    <- reticulate::py
  invisible(TRUE)
}

#' Union of probe_ids across all NNET checkpoints (used to extend imputation).
#'
#' @param nnet_files Character vector of .pth paths
#' @return Character vector of CpG IDs
get_nnet_required_cpgs <- function(nnet_files) {
  if (length(nnet_files) == 0) return(character(0))
  init_nnet_python()
  probes <- character(0)
  for (p in nnet_files) {
    info <- tryCatch(.nnet_py$py$load_nnet_checkpoint(normalizePath(p)),
                     error = function(e) NULL)
    if (!is.null(info)) {
      probes <- union(probes, as.character(info$probe_ids))
    }
  }
  probes
}

#' Run NNET inference for all checkpoints and return per-(sample, signature)
#' summarised scores using the same trim-min-trim-max-then-average rule as
#' the SVM pipeline (process_results).
#'
#' @param imputed_long Imputed beta data frame (long form with IlmnID column +
#'        one column per sample). Same shape returned by perform_imputation().
#' @param nnet_files Character vector of .pth paths
#' @param test_data_ids IDs of probands to keep
#' @return Tibble with columns SampleID, Signature, pNNET_average, pNNET_sd
make_nnet_predictions <- function(imputed_long, nnet_files, test_data_ids) {
  if (length(nnet_files) == 0) {
    return(tibble::tibble(SampleID = character(0),
                          Signature = character(0),
                          pNNET_average = numeric(0),
                          pNNET_sd = numeric(0),
                          n_nnet = integer(0)))
  }

  init_nnet_python()
  py <- .nnet_py$py

  # imputed_long: IlmnID + sample columns. Build a probe-indexed matrix.
  beta <- imputed_long
  rownames(beta) <- beta$IlmnID
  beta$IlmnID <- NULL
  beta <- beta[, intersect(test_data_ids, names(beta)), drop = FALSE]

  raw <- list()  # list of data.frames with (SampleID, Signature, Score, Key)

  for (mf in nnet_files) {
    info <- tryCatch(py$load_nnet_checkpoint(normalizePath(mf)),
                     error = function(e) {
                       warning("Failed to load NNET checkpoint ", basename(mf),
                               ": ", e$message)
                       NULL
                     })
    if (is.null(info)) next

    probes <- as.character(info$probe_ids)
    signature <- if (is.null(info$signature) || identical(info$signature, NA))
                   sub("^bayesian_([^_]+).*$", "\\1", tools::file_path_sans_ext(basename(mf)))
                 else as.character(info$signature)

    missing <- setdiff(probes, rownames(beta))
    if (length(missing) > 0) {
      warning(basename(mf), ": ", length(missing), "/", length(probes),
              " probes missing - skipping")
      next
    }

    sub_mat <- beta[probes, , drop = FALSE]
    keep_samples <- colnames(sub_mat)[colSums(is.na(sub_mat)) == 0]
    if (length(keep_samples) == 0) {
      warning(basename(mf), ": all samples have NaN in required probes - skipping")
      next
    }
    X <- t(as.matrix(sub_mat[, keep_samples, drop = FALSE]))  # n_samples x n_probes

    probs <- tryCatch(as.numeric(py$nnet_predict(info$model, X)),
                      error = function(e) {
                        warning(basename(mf), ": predict failed - ", e$message)
                        NULL
                      })
    if (is.null(probs)) next

    raw[[basename(mf)]] <- data.frame(
      SampleID  = keep_samples,
      Signature = signature,
      Score     = probs,
      Key       = basename(mf),
      stringsAsFactors = FALSE
    )
  }

  if (length(raw) == 0) {
    return(tibble::tibble(SampleID = character(0),
                          Signature = character(0),
                          pNNET_average = numeric(0),
                          pNNET_sd = numeric(0),
                          n_nnet = integer(0)))
  }

  out <- dplyr::bind_rows(raw)
  out <- out[out$SampleID %in% test_data_ids, , drop = FALSE]

  # Same trim rule as SVM process_results: drop the highest and lowest score
  # across replicate checkpoints, then mean/sd over the rest.
  # Capture raw checkpoint count per signature *before* the trim filter, so
  # the QC column reflects how many .pth checkpoints actually loaded for this
  # signature (independent of the trim rule that drops min+max).
  out %>%
    dplyr::group_by(SampleID, Signature) %>%
    dplyr::mutate(n_nnet = dplyr::n(),
                  rank   = rank(Score, ties.method = "first")) %>%
    dplyr::filter(n_nnet < 3 |
                   (rank != min(rank) & rank != max(rank))) %>%
    dplyr::summarise(pNNET_average = mean(Score),
                     pNNET_sd      = stats::sd(Score),
                     n_nnet        = dplyr::first(n_nnet),
                     .groups       = "drop")
}
