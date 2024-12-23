# Use the official Ubuntu base image
FROM ubuntu:20.04

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Update and install essential packages
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Create git clone
RUN git clone https://github.com/pedro-albino-duarte/VariGene.git /varigene-pedrod

# Install Miniconda
ENV PATH="/opt/conda/bin:$PATH"
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    conda update -n base -c defaults conda -y

# Add Conda to PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Configure Conda and install Bioconda
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda update -y conda

# Install Python (latest version)
RUN conda install -y python=3.11 && conda clean -a


# Download the reference genome hg19 (GRCh37) from UCSC
RUN wget -P /data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && \
    gunzip /data/hg19.fa.gz

# Consider the indexed files of the reference genomes
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.bwt"
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.sa"
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.ann"
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.pac"
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.0123"
RUN wget -P /data/ "https://storage.googleapis.com/bucket-challenge-pedrod/hg19.fa.amb"


# Set the working directory inside the container
WORKDIR /varigene-pedrod

# Install Numpy
RUN conda install -y numpy && conda clean -a

# Install Pandas
RUN conda install -y pandas && conda clean -a

# Install MultiQC
RUN conda install -y multiqc && conda clean -a

# Install Fastp
RUN conda install -y bioconda::fastp && conda clean -a

#Install BWA
RUN conda install -y bioconda::bwa && conda clean -a

#Install Seqkit
RUN conda install -y bioconda::seqkit && conda clean -a

#Install Samtools
RUN conda install -y bioconda::samtools && conda clean -a

#Install BCFtools
RUN conda install -y bioconda::bcftools && conda clean -a

# Set entry point
CMD ["tail", "-f", "/dev/null"]