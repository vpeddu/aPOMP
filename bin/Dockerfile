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
                   unzip && \
    pip3 install biopython \
                 ete3 \
		 cython \
		 pysam \
         typing \
         pickle5

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

