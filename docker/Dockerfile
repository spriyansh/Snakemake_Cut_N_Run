# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Metadata
LABEL base.image="ubuntu:20.04" \
    description="Cut & Run Data Analysis Custom Pipeline with Snakemake" \
    maintainer="<spriyansh29@gmail.com>"

# Avoid interactive Install
ENV DEBIAN_FRONTEND=noninteractive \
    R_VERSION=R-4.3.0 \
    PATH=$PATH:~/.local/bin/

# Update and Install necessary libraries, utilities and compilers
RUN apt-get update && \
    apt-get install -y \
    curl unzip wget gcc make openjdk-17-jdk xvfb \
    libxml2-dev libpng-dev gdal-bin libgdal-dev \
    libfontconfig1-dev libudunits2-dev libharfbuzz-dev \
    libfribidi-dev libglpk40 texlive texlive-latex-extra \
    texlive-fonts-recommended libxt6 \
    libcairo2-dev patch libgeos-dev libxt-dev \
    cmake libgsl-dev build-essential libfftw3-dev \
    libmagick++-dev libgsl-dev libssl-dev libxml2-dev \
    fastqc trim-galore bowtie2 samtools picard bedtools igv \
    git python3.8 python3-pip libcurl4-openssl-dev cmake \
    zlib1g-dev libbz2-dev liblzma-dev libxcb1 \
    libgsl-dev libssl-dev libxml2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --upgrade pip numpy MACS2 deeptools multiqc snakemake

# Download and setup bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig && \
    mv bedGraphToBigWig /usr/local/bin/

# Compile and install HTSLib
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 && \
    tar -jxvf htslib-1.14.tar.bz2 && \
    cd htslib-1.14 && \
    ./configure --prefix=/opt/htslib && \
    make && make install

# Compile and install preseq
WORKDIR /opt
RUN wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz && \
    tar -xzvf preseq-3.1.2.tar.gz
WORKDIR /opt/preseq-3.1.2
RUN ./configure LDFLAGS="-L/opt/htslib/lib" CPPFLAGS="-I/opt/htslib/include" --prefix=/opt/preseq --enable-hts && \
    make && make install

# Additional setups
WORKDIR /
COPY kseq_Compiled /additional_tools/
RUN ln -s /usr/bin/python3 /usr/local/bin/python && \
    git clone https://github.com/broadinstitute/picard.git && \
    cd picard && \
    ./gradlew shadowJar && \
    mkdir -p indexFolder rawFQ outputFolder intermediateFolder additional_tools

# Installing R
WORKDIR /usr/local
RUN wget https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz && \
    tar xvf $R_VERSION.tar.gz && \
    cd $R_VERSION && \
    ./configure --with-readline=no --with-x=no && \
    make -j 4 && \
    make install

# Clone, build, and install fastq_pair
WORKDIR /additional_tools/
RUN git clone https://github.com/linsalrob/fastq-pair.git && \
    cd fastq-pair && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# Download and install SEACR
WORKDIR /additional_tools/
RUN git clone https://github.com/FredHutch/SEACR.git

# Install BioConductor & DESeq2
RUN Rscript -e 'install.packages(c("BiocManager", "tidyverse", "ggpubr", "viridis", "corrplot"),repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'
RUN Rscript -e "BiocManager::install('DESeq2', update = TRUE, ask = FALSE, Ncpus = 24)"
RUN Rscript -e "BiocManager::install('chromVAR', update = TRUE, ask = FALSE, Ncpus = 24)"
RUN Rscript -e 'install.packages(c("ggcorrplot"),repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'


RUN Rscript -e "BiocManager::install(c('ChIPQC', 'ChIPseeker'), update = TRUE, ask = FALSE, Ncpus = 24)"
RUN Rscript -e "BiocManager::install(c('clusterProfiler', 'AnnotationDbi'), update = TRUE, ask = FALSE, Ncpus = 24)"
RUN Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'EnsDb.Hsapiens.v75'), update = TRUE, ask = FALSE, Ncpus = 24)"


WORKDIR /
RUN mkdir -p outputFolder indexFolder rawFQ scripts

WORKDIR /scripts

# Define ENTRYPOINT
ENTRYPOINT [ "/bin/bash" ]
