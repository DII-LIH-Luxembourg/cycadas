FROM rocker/shiny:latest

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    libgit2-dev \
    libhdf5-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    cmake \
    gfortran \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools', 'BiocManager', 'remotes'), repos='https://cloud.r-project.org/')"

RUN R -e "BiocManager::install(c('CATALYST', 'FlowSOM', 'SingleCellExperiment', 'ComplexHeatmap', 'flowCore', 'ConsensusClusterPlus'))"

RUN R -e "devtools::install_github('DII-LIH-Luxembourg/cycadas', dependencies = TRUE)"

EXPOSE 3838

CMD ["R", "-e", "options(shiny.host = '0.0.0.0', shiny.port = 3838); library(cycadas); cycadas()"]
