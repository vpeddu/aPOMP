FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3


RUN apt update && \
    apt install -y python3 && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
                   curl \
		           bedtools \
                   samtools \
                   pigz \
                   git \
                   bc \
                   unzip && \
    pip3 install biopython \
                 ete3 \
		 cython \
		 pysam \
         typing \
         pickle5 \
         taxopy \
         pyfasta

# minimap2 binary 
RUN wget https://github.com/lh3/minimap2/releases/download/v2.21/minimap2-2.21_x64-linux.tar.bz2 && \
    tar -xf minimap2-2.21_x64-linux.tar.bz2 && \
    mv minimap2-2.21_x64-linux/minimap2 /usr/bin/

RUN apt install -y zlib1g-dev 

RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure --prefix=/usr/bin/ --without-curses --disable-bz2 --disable-lzma  && \
    make && \
    make install 

# Install conda, bowtie2, and krakenuniq
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" && \
    conda install -c bioconda ucsc-fasplit seqkit blast hmmer libcurl seqtk

RUN git clone https://github.com/fenderglass/Flye && \
    cd Flye && \
    make

RUN mkdir amrfinder && cd amrfinder && \
# Get the URL of the latest binary tarball (linux only)
URL=`curl -s https://api.github.com/repos/ncbi/amr/releases/latest \
    | grep "browser_download_url.*amrfinder_binaries" \
    | cut -d '"' -f 4` && \
    # Download the latest binary tarball
curl -sOL "$URL" && \
# untar AMRFinderPlus
filename=`basename $URL` && \
tar xvfz $filename

# installing GNU parallel
RUN apt-get install -y parallel
RUN yes 'will cite' | parallel --citation
