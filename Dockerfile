# Use the official R Shiny base image
FROM rocker/shiny:latest

# Install system dependencies
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
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    cmake \
    gfortran \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install build tools
RUN R -e "install.packages(c('remotes', 'BiocManager', 'readr', 'zip'), repos='https://cloud.r-project.org/')"

# Install heavy Bioconductor dependencies first (this layer is cached, so it won't rebuild next time)
RUN R -e "BiocManager::install(c('CATALYST', 'FlowSOM', 'SingleCellExperiment', 'ComplexHeatmap', 'flowCore', 'ConsensusClusterPlus'))"
RUN apt-get update && apt-get install -y zip

# --- CHANGED SECTION STARTS HERE ---

# Copy the local project files into the container
COPY . /app

# Install the package from the local copy
# This works for private repos and avoids internet connection issues
RUN R -e "remotes::install_local('/app', dependencies = TRUE, upgrade = 'always')"

# Clean up source files to keep image small (optional)
RUN rm -rf /app

# --- CHANGED SECTION ENDS HERE ---

# Expose the port
EXPOSE 3838

# Start the app
CMD ["R", "-e", "options(shiny.host = '0.0.0.0', shiny.port = 3838); library(cycadas); cycadas()"]
