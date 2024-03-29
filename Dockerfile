FROM rocker/r-ver:4.1.2

MAINTAINER Alexandre Schickele "alexandre.schickele@imev-mer.fr"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    libudunits2-dev \
    libproj-dev \
    libgeos-dev \
    libgdal-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git 

# zlib for R package httpuv
RUN apt-get -y install zlib1g-dev

# R CRAN packages
RUN Rscript -e "install.packages(c('devtools'), repos='https://cran.r-project.org/')"
RUN Rscript -e "install.packages(c('XML'), repos='https://cran.r-project.org/')"
RUN Rscript -e "install.packages(c('shiny','shinyjs','shinythemes'), repos='https://cran.r-project.org/')"

RUN Rscript -e "install.packages(c('raster','virtualspecies'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN Rscript -e "install.packages(c('RColorBrewer','colorspace'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN Rscript -e "install.packages(c('tidyverse'), repos='https://cran.r-project.org/', dependencies=TRUE)"

# Copy application code and data
COPY app .

# Define command to run docker at start
CMD ["R", "-e shiny::runApp('.',port=3838,host='0.0.0.0')"]

# END
EXPOSE 3838