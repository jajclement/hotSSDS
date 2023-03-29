FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="pauline.auffret@ifremer.fr" \
      description="Docker image containing peak calling tools requirements for SSDS DMC1 pipeline"

COPY environment.yml /
RUN micromamba create -f /environment.yml
RUN micromamba clean -a
ENV PATH /opt/conda/envs/peak-calling-box/bin:$PATH

