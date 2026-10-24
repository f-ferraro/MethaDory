FROM ghcr.io/prefix-dev/pixi:noble AS build

WORKDIR /app
COPY pixi.toml .
COPY pixi.lock .
RUN pixi install --locked --environment default
RUN echo "#!/bin/bash" > /app/entrypoint.sh && \
    pixi shell-hook --environment default -s bash >> /app/entrypoint.sh && \
    echo 'exec "$@"' >> /app/entrypoint.sh

FROM ghcr.io/prefix-dev/pixi:noble AS final

WORKDIR /app
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /app/pixi.toml /app/pixi.toml
COPY --from=build /app/pixi.lock /app/pixi.lock
# The ignore files are needed for 'pixi run' to work in the container
COPY --from=build /app/.pixi/.gitignore /app/.pixi/.gitignore
COPY --from=build /app/.pixi/.condapackageignore /app/.pixi/.condapackageignore
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh
COPY ./data /app/data
COPY ./src /app/src
COPY ./ONT-Input-preparation /app/ONT-Input-preparation

RUN echo "#!/bin/bash" > /app/entrypoint.sh && \
    echo "pixi run \"\$@\"" >> /app/entrypoint.sh

# run it to install the missing dependencies.
RUN pixi run Rscript -e "BiocManager::install('ChAMPdata')"
RUN pixi run Rscript -e "BiocManager::install('FDb.InfiniumMethylation.hg19')"
RUN pixi run Rscript -e "BiocManager::install('GenomeInfoDb')"
RUN pixi run Rscript -e "BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19')"
# hexylena: unknown why this was needed a second time, but failed otherwise during --help to install ChAMPdata
RUN pixi run Rscript -e "BiocManager::install('ChAMPdata')"
# hexylena: and still required, otherwise further installation commands run on execution.
RUN pixi run MethaDory_cli --help

ENTRYPOINT [ "/app/entrypoint.sh" ]
