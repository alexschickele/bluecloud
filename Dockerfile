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

#R GitHub packages (with release in CRAN, but not re-published yet)

#R GitHub packages (not yet released in CRAN)
 
EXPOSE 3838

RUN apt-get install -y curl
CMD ["R", "-e shiny::runApp('/srv/shiny/bluecloud_wp3_d2',port=3838,host='0.0.0.0')"]
