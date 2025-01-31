# pull base image
FROM rocker/tidyverse:4.2.2

# docker / podman build -t cancerbits/dockr:noorizadeh_rstudio - < noorizadeh_rstudio.Dockerfile

LABEL maintainer Florian Halbritter "florian.halbritter@ccri.at"
LABEL version n-4.2.2

# change some permissions
RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

RUN apt-get -y update && apt-get -y install \
  gcc \
  libncurses5-dev \
  libncursesw5-dev \
  liblzma-dev \
  bzip2 \
  libbz2-dev \ 
  libglpk40 \
  libxt6 \
  libgsl-dev \
  imagemagick \
  libmagick++-dev \
  libmagickcore-dev \
  && apt-get clean \
  && rm -rf /tmp/* /var/tmp/*
 
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \	
	tar xvfj samtools-1.16.1.tar.bz2 && \ 
	cd samtools-1.16.1 && \
	./configure && \
	make && \
	make install

# The default CRAN mirror and Bioconductor versions
ENV CRAN=https://packagemanager.posit.co/cran/2022-04-21
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2022-04-21'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site
RUN R -e "BiocManager::install(version = "3.16")"

# Install R libraries
RUN R -e "BiocManager::install(c('Rhtslib', 'rtracklayer', 'Rsamtools', 'BSgenome', 'GenomicFeatures', 'Rsubread'))"
RUN R -e "BiocManager::install(c('markdown', 'patchwork', 'R.utils', 'viridis', 'ggVennDiagram', 'eulerr', 'dendsort', 'ggplot2', 'data.table', 'ggrepel'))"
RUN R -e "BiocManager::install(c('simpleCache', 'ComplexHeatmap', 'magick', 'ggrastr', 'fastcluster', 'MKmisc', 'ggplotify', 'stringr'))"


RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'JASPAR2020', 'motifmatchr', 'TFBSTools'))"
RUN R -e "BiocManager::install(c('hypeR','fgsea'))"


RUN R -e "remotes::install_github(repo = 'cancerbits/canceRbits@v0.1.6', upgrade='never')"

RUN R -e "devtools::install_github('da-bar/JASPAR2022@5e146566d12646e275e51d210d87ae328f9c244d')"


RUN R -e "BiocManager::install(c('org.Mm.eg.db', 'EnhancedVolcano', 'ggvenn', 'msigdbr'))" 

RUN apt-get -y update && apt-get -y install \
  ncftp \
  && apt-get clean \
  && rm -rf /tmp/* /var/tmp/*
