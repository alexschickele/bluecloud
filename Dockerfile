FROM rocker/r-ver:4.0.5

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
RUN R -e "install.packages(c('shiny','shinydjs','shinythemes'), repos='https://cran.r-project.org/')"

RUN R -e "install.packages(c('raster'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('tidyverse'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('leaflet'), repos='https://cran.r-project.org/', dependencies=TRUE)"

# COPY DATA
WORKDIR home/aschickele/workspace/bluecloud
COPY app .

# END
EXPOSE 3838

RUN apt-get install -y curl
CMD ["R", "-e shiny::runApp('/srv/shiny/bluecloud_wp3_d2_shiny',port=3838,host='0.0.0.0')"]
