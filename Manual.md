# MethaDory manual

This manual explains how to read the results and figures produced by MethaDory.

---

## 1. Outputs

### Tables

The `.xlsx` workbook written by `MethaDory_html` contains:

| Sheet | Content |
|---|---|
| `Analysis_Info` | Signature version and file used, and when the run was made, so results can be traced to the exact signature set |
| `Predictions_Long` | One row per sample and signature (columns below) |
| `Predictions_Wide` | The same values with one row per sample and the signatures spread across columns |
| `Signatures_Summary` | Per signature, the highest, lowest and mean `pCombined` across the samples of the run, and how many samples fall in each confidence bin (section 2.1) |
| `Cell_Proportions` | Estimated blood cell fractions per sample (when sample statistics are computed) |
| `Methylation_Age` | Predicted age and the number of clock CpGs that were missing |
| `Chromosomal_Sex` | X and Y signal and the predicted sex |

The app's *Download Results* workbook holds the same numbers in a simpler layout, for the probands selected in the sidebar: `Analysis Info`, `MethaDory Predictions` (wide: one row per sample, with the combined score under its internal name `mean case`), `Cell deconvolutions`, `Predicted age` and `Predicted chr sex`. 

Columns of the prediction table:

| Column | Meaning |
|---|---|
| `SampleID` | Sample name, from the input header |
| `SVM` | Signature name |
| `pSVM_average`, `pSVM_sd` | Mean and standard deviation of the SVM probability across the SVM models trained for this signature |
| `pNNET_average`, `pNNET_sd` | The same for the neural networks |
| `pCombined` | Mean of `pSVM_average` and `pNNET_average`. This is the score used for ranking, thresholds and confidence bins. If a signature has models from only one family, the other average is empty and `pCombined` equals the available one; the prediction plot then shows no whisker |
| `n_svm`, `n_nnet` | Number of models that contributed to each average |
| `pct_na_pre` | Percentage of this signature's probes that were missing in the input before imputation |



---

### The prediction plot

- One position on the x-axis per signature; the y-axis is the score from 0 to 1.
- The **point** is `pCombined`.
- The **whiskers** run from the lower to the higher of the two family scores (`pSVM_average`, `pNNET_average`). They show disagreement between the two classifiers, not a confidence interval.
- Grey guide lines mark 0, 0.5 and 1. A red dashed line marks the display threshold when one is set.
- When several samples are selected, each is drawn in a different colour.

How to read it:

| Pattern | Reading |
|---|---|
| High point, short whisker | Both classifier families agree that the signature is present |
| High point, long whisker | One family scores high and the other low. Treat as unresolved and rely on the per-signature figure |
| Several related signatures elevated together | Expected for signatures that share probes or biology (for example genes of the same complex). Compare their per-signature figures |
| Many unrelated signatures elevated | Could be caused by a sample problem  (check the QC panels and `pct_na_pre`), newborns have different cell composition, samples that are not blood, very strong DNAm alterations due to other disoders or drug treatments |
| Everything near 0 | No evidence for any tested signature. This does not exclude a disorder, the variant may not produce the signature, or the signature may not be in the collection |

---

## Interpreting the sample QC panels

These panels tell you whether the sample can be trusted at all. Look at them before the scores and the per-signature figures.

### Missing values before imputation

CpGs absent from the input file count as missing and are imputed. The more values are imputed, the more the prediction depends on the imputation reference rather than on the sample. The per-signature equivalent is `pct_na_pre` in the prediction table.

### Methylation age

The table gives the predicted age and the number of clock CpGs that were missing. Compare it with the recorded age as an identity and quality check. A large unexplained difference is a reason to check the sample. Some disorders are associated with altered methylation age, so a difference is not automatically an error.

### 3.3 Blood cell proportions

- **Plot.** Grey points are the estimated cell fractions of the samples the classifiers were trained on; large coloured points are the selected samples.
- **Table.** Each cell type of each sample is compared with the training distribution:

| Status | Meaning |
|---|---|
| `PASS` | Within 3 standard deviations of the training mean |
| `WARNING` | More than 3 standard deviations from the mean, but still inside the observed training range |
| `FAIL` | Outside the range observed in training |


Episignatures are measured against a background of normal blood composition. A sample with unusual composition (infection, treatment, haematological disease, a tissue other than blood) will have unexpected results.

### 3.4 Chromosomal sex

The plot places each sample by its X-chromosome signal (x-axis) and Y-chromosome signal (y-axis). 
Compare the predicted sex with the recorded sex. A mismatch could suggest a sample swap, contamination, or a sex chromosome aneuploidy. ONT and PacBio XY samples are sometimes predicted as XXY.

### 3.5 PCA against controls

Principal component analysis of each selected sample together with the control samples used for developing the classifiers, computed on the beta values before imputation. Only autosomal CpGs measured in the sample are used, and of those the top 1% most variable. Two panels are shown per sample: PC1 vs PC2 and PC3 vs PC4. Controls are grey and the sample is red.

A sample lying far from the controls may be of poor quality, or come from a different tissue or platform, and its predictions should be interpreted with caution.

### 3.6 Beta-value distribution against controls

Distribution of the beta values before imputation over all autosomal CpGs measured in the sample. Controls are grey and the sample is red.

Controls show two peaks, near 0 and near 1. A distribution that departs from this shape may indicate poor quality, or data from a platform that is not on the array scale (for example uncorrected long-read calls).

---

## 4. Interpreting the per-signature figure

One figure is drawn for every signature at or above the score threshold (section 2.1). It uses only that signature's CpGs and shows the proband together with a reference set.

### 4.1 Who is in the figure

| Group | Colour | What it is |
|---|---|---|
| Proband | pink, drawn larger or as a triangle | The selected input sample(s) |
| Real cases | orange, labelled with the diagnosis | Affected individuals from `data/affectedindividuals` |
| In-silico cases | pale yellow | Synthetic cases: a control profile plus the signature's effect. Added only when there are fewer than N real cases |
| Controls | teal | The N closest controls |

### 4.2 Layout

| Row | Left | Right |
|---|---|---|
| 1 | PCA | PCA with the proband projected |
| 2 | MDS | Delta concordance |
| 3 | Similarity to median profiles | Sample–sample correlation heatmap |
| 4 | Ranked neighbours  | |
| 5 | Probe heatmap | |

### PCA and PCA with the proband projected

**Left panel** Principal component analysis of all displayed samples, proband included. Axis labels give the percentage of variance explained. The proband falling inside the case cloud supports the score

**Right panel** Because the proband takes part in the fit in the PCA on the left, a proband from a different platform can end up far from both groups. 
In the right panel, the components are computed from controls and cases only, and the proband is then projected onto them.  If the proband moves a lot between the left and right panels, the left panel was being driven by the proband's own noise.

### 4.5 MDS (row 2, left)

An ordination like PCA, but built from sample-to-sample distances (1 − Spearman rank correlation) computed on the probes each pair of samples shares, after subtracting the control median from every probe.

### 4.6 Delta concordance (row 2, right)

One point per CpG.

- **x-axis:** the signature's effect at that CpG (median case − median control).
- **y-axis:** the proband's deviation at that CpG (proband − median control).
- **Orange diagonal:** where a typical case lies (slope 1). **Teal horizontal:** where a control lies (slope 0).
- **Pink line:** the fit through the proband's points. The subtitle gives its slope, offset and correlation `r`.

| Quantity | Reading |
|---|---|
| Slope near 1 | The proband carries the full signature |
| Slope near 0 | The proband resembles a control |
| Slope in between | Attenuated signature: mosaicism, a variant with a milder effect, or a related disorder that shares part of the signature |
| Offset away from 0 | A global shift that is the same at every CpG, typically platform or batch. It is separate from the signature and does not count for or against it |
| Low `r` | The points do not follow a line: whatever the slope, the evidence is weak |

### 4.7 Similarity to median profiles (row 3, left)

Pearson correlation of each sample with the median control profile (x) and with the median case profile (y). Above the red diagonal means more similar to cases; below means more similar to controls.

### 4.8 Sample–sample correlation heatmap (row 3, right)

Spearman (rank) correlation between every pair of displayed samples, clustered. Before correlating, the median of the displayed controls is subtracted from every probe, so each sample is described by its **deviation from the control profile** rather than by its raw beta values. Raw values share each probe's baseline methylation level, which would make every pair of samples correlate close to 1 and hide the groups. The top bar gives each sample's group; the proband is labelled on the right.

### 4.9 Ranked neighbours (row 4)

The proband, then all displayed reference samples ordered by their distance to it, nearest on the left. The proband itself sits at rank 0 and distance 0 and is the sample every other one is measured from, so the height of the first reference point reads directly as "how far is the nearest sample". The colour bar gives each sample's group; the points below give the actual distances. The panel title counts how many of the N nearest samples are cases, where N is the number of cases shown.

- Case colours stacked on the left with a jump in distance before the first control: the proband's neighbourhood is cases.
- Colours mixed: no clear neighbourhood.

Unlike a dendrogram, this view is centred on the proband, has no arbitrary ordering, and shows how large the differences are.

### 4.10 Probe heatmap (row 5)

Rows are the signature's CpGs, columns are samples; both are clustered. Values are standardised per CpG: yellow is above that CpG's mean, blue is below. The annotation bars give group, platform, sex and age group; the proband is labelled underneath.

- Look at the proband's column: does it show the same blocks of yellow and blue as the cases, across the whole signature or only over part of it?
- Look at the annotation bars: if the samples group by platform, sex or age rather than by case status, the clustering is driven by a confounder.

---

## 5. Putting it together

1. **Sample QC first.** A high share of missing values, a failed cell composition check, a sex mismatch, an implausible age, or a sample lying away from the controls in the QC PCA or beta distribution all lower confidence in every score for that sample.
2. **Scores.** Note which signatures are elevated, whether the two classifier families agree, and how many probes had to be imputed (`pct_na_pre`).
3. **Per-signature figure.** For each elevated signature:
   - Is the proband with the cases in the projected PCA and the MDS?
   - Is the delta-concordance slope close to 1, with a good `r`?
   - Are its nearest neighbours cases, ideally real ones?
   - Does the probe heatmap show the case pattern across the signature?
4. **Weigh agreement.** A signature is best supported when a clean sample, a high score with agreeing classifiers, and consistent per-signature panels all point the same way. When the panels disagree with each other or with the score, treat the result as inconclusive.

---

## 6. Caveats

- **Missing data.** A high proportion of missing data increases run time and can reduce classifier performance. The percentage of signature probes that were missing before imputation is reported for every sample and signature (`pct_na_pre`).
- **Array normalisation.** The array classifiers were trained on samples normalised with Noob. Performance varies with other normalisation methods, so please use Noob.
- **Sample type.** The classifiers and the reference cohort are based on blood. Other tissues are out of scope and their results should be interpreted cautiously.
- **MethaDory provides supporting evidence.** Results should be interpreted by trained individuals together with the clinical picture and genetic findings, and are not a standalone diagnostic test.
- **A negative result does not exclude a disorder.** Not every pathogenic variant produces the signature, and only the signatures in the collection are tested.
- **In-silico cases are synthetic.** Each is a control profile plus the signature effect, so a group of them is more uniform than real patients would be. Panels that depend on the spread of the cases are optimistic when few real cases are available.
- Some signatures do not replicate in MethaDory, the current list is CHD8 ArefEsghi 2020, PURA Xiao 2024, Cul3 van der Laan 2025, KDM4B Levy 2022, CHD3 Santini 2025.
- **Overlapping and undiscovered signatures.** Related disorders can share probes and biology. Several elevated scores may describe one underlying pattern rather than several diagnoses.
- **MethaDory is provided as is, without warranty of any kind, and the authors accept no responsibility for decisions made on the basis of its results (see `LICENSE`).**


