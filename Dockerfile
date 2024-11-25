FROM rocker/shiny:4.3.2

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev libglpk-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN mkdir /home/data \
    && chown -R shiny:shiny /home/data
COPY /data/ /home/data/

COPY /start-script.sh ./start-script.sh
RUN chmod +x ./start-script.sh

RUN rm -rf /srv/shiny-server/*
COPY /app/ /srv/shiny-server/

# Command to install standard R packages from CRAN; enter the list of required packages for your app here
RUN Rscript -e 'install.packages(c("shiny","tidyverse","plotly","BiocManager", "Seurat", "DT", "dplyr", "ggplot2", "SeuratObject"), dependencies=TRUE)'

# Command to install packages from Bioconductor; enter the list of required Bioconductor packages for your app here
RUN Rscript -e 'BiocManager::install(c("Biostrings"),ask = F)'

USER shiny

EXPOSE 3838

ENTRYPOINT ["./start-script.sh"]
