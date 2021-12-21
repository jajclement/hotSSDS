FROM nfcore/base:1.9
MAINTAINER Pauline Auffret <pauline.auffret@igh.cnrs.fr>
LABEL authors="Pauline Auffret" \
    description="Docker image containing all requirements for SSDS nextflow pipeline"

# Install conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/ssdsnextflowpipeline/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ssdsnextflowpipeline_v2.0 > ssdsnextflowpipeline_v2.0.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

