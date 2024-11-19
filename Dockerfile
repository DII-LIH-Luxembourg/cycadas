# Use an R base image or a Linux image
FROM rocker/shiny

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    && apt-get clean

# Install R dependencies
RUN R -e "install.packages('devtools')" \
    && R -e "install.packages(c('shiny', 'dplyr', 'ggplot2', 'BiocManager'))" \
    && R -e "BiocManager::install(c('flowCore', 'mousetrap'))"


# Set the working directory
WORKDIR /usr/src/app
COPY . /usr/src/app

RUN R -e "devtools::install('/usr/src/app')"
# Expose the Shiny server port
EXPOSE 3838

# Start the CyCadas Shiny app
CMD ["R", "-e", "cycadas::cycadas()"]

