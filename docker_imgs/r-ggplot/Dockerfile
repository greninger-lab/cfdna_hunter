FROM python:3.11

RUN pip3 install plotly==5.18.0 
FROM r-base

# install packages
# these are ones I like
RUN echo 'install.packages(c("ggplot2", "plyr", "reshape2", "RColorBrewer", "scales","grid", "wesanderson", "svglite"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
    && Rscript /tmp/packages.R

RUN apt-get update
RUN apt-get install -y procps

