FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="pauline.auffret@ifremer.fr" \
      description="Docker image containing bigwig manipulation tools requirements for SSDS DMC1 pipeline"

COPY environment.yml /
RUN micromamba create -f /environment.yml
RUN micromamba clean -a
ENV PATH /opt/conda/envs/bigwig-box/bin:$PATH

