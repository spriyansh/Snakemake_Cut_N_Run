# Use rocker/rstudio as the base image
# Build: docker build --file Rstudio.Dockerfile --cpu 24 -t spriyansh29/cnr_snake_rstudio .
FROM rocker/rstudio:4.3

# Maintainer
LABEL maintainer="Priyansh Srivastava <spriyansh29@gmail.com>"

# Install Required OS libs
RUN apt-get update && \
    apt-get install -y \
        libxml2-dev \
        libpng-dev \
        gdal-bin \
        libgdal-dev \
        libfontconfig1-dev \
        libudunits2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libglpk40 \
        texlive \
        texlive-latex-extra \
        texlive-fonts-recommended \
        libxt6 \
        libcairo2-dev \
        patch \
        libgeos-dev \
        libxt-dev \
        cmake \
        libgsl-dev \
        build-essential \
        libfftw3-dev \
        libmagick++-dev \
        libgsl-dev \
        libssl-dev \
        libxml2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R Packages and BioConductor libraries
RUN Rscript -e 'install.packages(c("tidyverse", "BiocManager"), repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)' && \
    Rscript -e "BiocManager::install(c('DESeq2', 'chromVAR'), update = TRUE, ask = FALSE, Ncpus = 24)"

RUN Rscript -e 'install.packages(c("ggpubr", "styler", "viridis"), repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'

RUN Rscript -e 'install.packages(c("corrplot"), repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'

RUN Rscript -e "BiocManager::install(c('ChIPQC', 'ChIPseeker'), update = TRUE, ask = FALSE, Ncpus = 24)"

RUN Rscript -e "BiocManager::install(c('clusterProfiler', 'AnnotationDbi'), update = TRUE, ask = FALSE, Ncpus = 24)"

RUN Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'EnsDb.Hsapiens.v75'), update = TRUE, ask = FALSE, Ncpus = 24)"

# Create ESR data directory
RUN mkdir /home/priyansh
RUN mkdir /home/Rscripts

# Expose port 8787 (RStudio runs on this port)
EXPOSE 8787

CMD ["/init"]
