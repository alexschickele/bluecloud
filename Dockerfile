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
   
RUN apt-get update && apt-get upgrade -y

# R CRAN packages
RUN R -e "install.packages(c('devtools'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('XML'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('shiny','shinyjs','shinythemes'), repos='https://cran.r-project.org/')"

RUN R -e "install.packages(c('raster','virtualspecies'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('RColorBrewer','colorspace'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('tidyverse'), repos='https://cran.r-project.org/', dependencies=TRUE)"

# COPY DATA
COPY app .

# END
EXPOSE 3838

RUN apt-get install -y curl
CMD ["R", "-e shiny::runApp('.',port=3838,host='0.0.0.0')"]
