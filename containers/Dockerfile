# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:${PATH}"

# Update package list and install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libgsl-dev \
    libgslcblas0 \
    libhdf5-dev \
    libboost-all-dev \
    cmake \
    pkg-config \
    unzip \
    default-jre \
    default-jdk \
    r-base \
    r-base-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-setuptools \
    python3-wheel \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages (only scientific packages)
RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    pysam \
    pyvcf \
    umi_tools

# Install Java (required for GATK)
RUN apt-get update && apt-get install -y openjdk-11-jdk

# Install BWA-MEM2
RUN cd /tmp && \
    wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    cp bwa-mem2-2.2.1_x64-linux/bwa-mem2 /usr/local/bin/ && \
    cd / && rm -rf /tmp/bwa-mem2-2.2.1_x64-linux*

# Install Samtools
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -xjf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure && make && make install && \
    cd / && rm -rf /tmp/samtools-1.16.1*

# Install BEDtools
RUN cd /tmp && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -xzf bedtools-2.30.0.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /usr/local/bin/ && \
    cd / && rm -rf /tmp/bedtools2*

# Install FastQC
RUN cd /tmp && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    cp FastQC/fastqc /usr/local/bin/ && \
    cd / && rm -rf /tmp/FastQC*

# Install Fastp
RUN cd /tmp && \
    wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.tar.gz && \
    tar -xzf v0.23.4.tar.gz && \
    cd fastp-0.23.4 && \
    make && \
    cp fastp /usr/local/bin/ && \
    cd / && rm -rf /tmp/fastp-0.23.4*

# Install GATK
RUN cd /tmp && \
    wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && \
    unzip gatk-4.2.6.1.zip && \
    cp gatk-4.2.6.1/gatk /usr/local/bin/ && \
    cd / && rm -rf /tmp/gatk-4.2.6.1*

# Install Tabix
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/1.16.1/htslib-1.16.1.tar.bz2 && \
    tar -xjf htslib-1.16.1.tar.bz2 && \
    cd htslib-1.16.1 && \
    ./configure && make && make install && \
    cd / && rm -rf /tmp/htslib-1.16.1*

# Install BCFtools
RUN cd /tmp && \
    wget https://github.com/samtools/bcftools/releases/download/1.16.1/bcftools-1.16.1.tar.bz2 && \
    tar -xjf bcftools-1.16.1.tar.bz2 && \
    cd bcftools-1.16.1 && \
    ./configure && make && make install && \
    cd / && rm -rf /tmp/bcftools-1.16.1*

# Install FGBio (for UMI consensus)
RUN cd /tmp && \
    wget https://github.com/fulcrumgenomics/fgbio/releases/download/2.0.2/fgbio-2.0.2.jar && \
    mv fgbio-2.0.2.jar /usr/local/bin/fgbio.jar && \
    echo '#!/bin/bash\njava -jar /usr/local/bin/fgbio.jar "$@"' > /usr/local/bin/fgbio && \
    chmod +x /usr/local/bin/fgbio && \
    cd / && rm -rf /tmp/*

# Install Mosdepth (for coverage analysis)
RUN cd /tmp && \
    wget https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth && \
    chmod +x mosdepth && \
    mv mosdepth /usr/local/bin/ && \
    cd / && rm -rf /tmp/*

# Install Bismark (for methylation analysis)
RUN cd /tmp && \
    wget https://github.com/FelixKrueger/Bismark/archive/0.24.1.tar.gz && \
    tar -xzf 0.24.1.tar.gz && \
    cd Bismark-0.24.1 && \
    cp bismark* /usr/local/bin/ && \
    cp deduplicate_bismark /usr/local/bin/ && \
    cp bismark_methylation_extractor /usr/local/bin/ && \
    cd / && rm -rf /tmp/Bismark-0.24.1*

# Install Bowtie2 (required by Bismark)
RUN cd /tmp && \
    wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-linux-x86_64.zip && \
    unzip bowtie2-2.5.1-linux-x86_64.zip && \
    cp bowtie2-2.5.1-linux-x86_64/bowtie2* /usr/local/bin/ && \
    cd / && rm -rf /tmp/bowtie2-2.5.1-linux-x86_64*

# Install ichorCNA dependencies
RUN R -e "install.packages(c('optparse', 'ggplot2', 'gplots', 'copynumber'), repos='https://cran.rstudio.com/')"

# Install ichorCNA
RUN cd /tmp && \
    git clone https://github.com/broadinstitute/ichorCNA.git && \
    cp -r ichorCNA /usr/local/ && \
    echo '#!/bin/bash\nRscript /usr/local/ichorCNA/scripts/runIchorCNA.R "$@"' > /usr/local/bin/runIchorCNA.R && \
    chmod +x /usr/local/bin/runIchorCNA.R && \
    cd / && rm -rf /tmp/ichorCNA

# Install Nextflow
RUN cd /tmp && \
    curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    cd / && rm -rf /tmp/*

# Set working directory
WORKDIR /workspace

# Default command
CMD ["/bin/bash"]
