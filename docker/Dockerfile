FROM rocker/r-ver:4.3.1

LABEL maintainer "Diego Garrido-Martín <dgarrido@ub.edu>" \
      version "0.0.1" \
      description "Docker image for digenic-nf"

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# bcftools dependencies
       libbz2-dev \
       liblzma-dev \
# python3
       python3 \
       python3-pip \
# nextflow trace (ps)
       procps \
# other
       wget \
       bc \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Set default CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >>"${R_HOME}/etc/Rprofile.site"

## Install bcftools and HTSlib
RUN wget --quiet https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 \
    && tar jxf bcftools-1.18.tar.bz2 \
    && cd bcftools-1.18 \
    && ./configure \
    && make \
    && make install

RUN cd bcftools-1.18/htslib-1.18 \
    && ./configure \
    && make \
    && make install

## Install PLINK2
RUN wget --quiet https://s3.amazonaws.com/plink2-assets/alpha4/plink2_linux_x86_64_20230922.zip \
    && unzip plink2_linux_x86_64_20230922.zip \
    && mv plink2 /usr/local/bin/

## Install python packages
RUN pip3 install numpy

## Install R packages
RUN R -e 'install.packages(c("optparse", "data.table"))'
