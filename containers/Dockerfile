FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
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
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libopenblas-dev \
    liblapack-dev \
    gfortran \
    pkg-config \
    cmake \
    unzip \
    bzip2 \
    gzip \
    tar \
    default-jre \
    default-jdk \
    && rm -rf /var/lib/apt/lists/*

# Install Python 3.8
RUN apt-get update && apt-get install -y python3.8 python3.8-dev python3-pip

# Set Python 3.8 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1

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

# Set working directory
WORKDIR /data

# Create directories for pipeline
RUN mkdir -p /data/{refs,resources,wes,plasma,pon,features,results,reports}

# Set default command
CMD ["/bin/bash"]
