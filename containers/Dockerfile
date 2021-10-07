## base container
FROM mambaorg/micromamba:0.15.3 as base_container
USER root
RUN apt-get update \
    && apt-get install --yes rename tini procps curl \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER micromamba
ENTRYPOINT ["tini", "--"]
CMD ["/bin/bash"]

## main container
FROM base_container

# custon conda file. Default: env.base.yml
ARG CONDA_FILE=containers/env.base.yml

# adding opencontainer labels to link registry to github repository
LABEL org.opencontainers.image.title="quick-rnaseq-nf"
LABEL org.opencontainers.image.description="A basic and quick workflow for differential expression analysis"
LABEL org.opencontainers.image.url="https://github.com/stracquadaniolab/quick-rnaseq-nf"
LABEL org.opencontainers.image.documentation="https://github.com/stracquadaniolab/quick-rnaseq-nf"
LABEL org.opencontainers.image.source="https://github.com/stracquadaniolab/quick-rnaseq-nf"
LABEL org.opencontainers.image.vendor="stracquadaniolab"
LABEL org.opencontainers.image.authors="Giovanni Stracquadanio"
LABEL org.opencontainers.image.revision="v0.2.0"

COPY --chown=micromamba:micromamba ${CONDA_FILE} /tmp
RUN micromamba install -y -n base -f /tmp/`basename ${CONDA_FILE}` && \
    micromamba clean --all --yes
