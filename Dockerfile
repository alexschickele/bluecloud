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

#Python Github libraries
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpq-dev python3.8 python3-pip python3-setuptools python3-dev

RUN pip install --upgrade pip
RUN pip install -e vcs+protocol://https://github.com/alexschickele/bluecloud/function/MBTR/#egg=pkg&subdirectory=pkg_dir

# install dependencies of the app

#R CRAN packages
RUN R -e "install.packages(c('devtools'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('XML'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('shiny'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('shinyjs'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('shinydashboard'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('jsonlite'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('DT'), repos='https://cran.r-project.org/')"
RUN R -e "install.packages(c('tibble'), repos='https://cran.r-project.org/')"

RUN R -e "install.packages(c('raster'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('abind'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('feather'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('reticulate'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('RColorBrewer'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('parallel'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('mvrsquared'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('tidyverse'), repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages(c('RSQLite'), repos='https://cran.r-project.org/', dependencies=TRUE)"

#R GitHub packages (with release in CRAN, but not re-published yet)
RUN R -e "devtools::install_github('eblondel/ows4R')"

EXPOSE 3838

RUN apt-get install -y curl
CMD ["R", "-e shiny::runApp('/srv/shiny/bluecloud_wp3_d2_shiny',port=3838,host='0.0.0.0')"]
