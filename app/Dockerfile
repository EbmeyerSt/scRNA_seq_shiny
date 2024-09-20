FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev libglpk-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Command to install standard R packages from CRAN; enter the list of required packages for your app here
RUN Rscript -e 'install.packages(c("shiny","tidyverse","BiocManager", "Seurat", "DT", "dplyr", "ggplot2", "SeuratObject"), dependencies=TRUE)'

# Command to install packages from Bioconductor; enter the list of required Bioconductor packages for your app here
RUN Rscript -e 'BiocManager::install(c("Biostrings"),ask = F)'

RUN rm -rf /srv/shiny-server/*
COPY ./ /srv/shiny-server/

# Make sure permissions are set so that the shiny user can read/write to the app directory.
# This assumes that there is a 'data' directory within './' that you reference in your Shiny app.
RUN chown -R shiny:shiny /srv/shiny-server && \
    chmod -R 755 /srv/shiny-server

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
