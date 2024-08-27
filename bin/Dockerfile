FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3


RUN apt update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-add-repository universe && \
    apt update && \
    apt-get clean

RUN   apt install -y wget \
                   curl \
                   bedtools \
                   pigz \
                   git \
                   bc \
                   unzip \
                   build-essential \
                   checkinstall \
                   libreadline-gplv2-dev \
                   libncursesw5-dev \
                   libssl-dev \
                   libsqlite3-dev \
                   tk-dev \
                   libgdbm-dev \
                   libc6-dev \
                   libbz2-dev \
                   liblzma-dev \
                   libffi-dev

RUN wget https://www.python.org/ftp/python/3.12.4/Python-3.12.4.tgz && \
    tar -xzf Python-3.12.4.tgz && \
    cd Python-3.12.4 && \
    ./configure && \
    make install

    #&& \
    #make altinstall

RUN apt install -y python-biopython \
                   python3-pip \
                   python3-pysam 

# Register the version in alternatives
#RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.12 1

# Set python 3 as the default python
#RUN update-alternatives --set python /usr/bin/python3.12

# Upgrade pip to latest version
# RUN curl -s https://bootstrap.pypa.io/pip/3.12/get-pip.py -o get-pip.py && \
#     python get-pip.py --force-reinstall && \
#     rm get-pip.py


RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" 

RUN conda install -c bioconda -c conda-forge \
    #ucsc-fasplit \
    seqkit \
    blast \
    hmmer\
    seqtk \
    krakenuniq 

# RUN /usr/local/miniconda/bin/python3 -m pip install cython
RUN /usr/local/miniconda/bin/python3  -m  pip install biopython \
         ete3 \
         typing \
         pickle5 \
         taxopy \
         pyfasta \
         ncbi-taxonomist \
         ete3 \
         pandas

RUN python3 -m pip install cython
RUN python3  -m  pip install biopython \
         ete3 \
         typing \
         #pickle5 \
         taxopy \
         pyfasta \
         ncbi-taxonomist \
         ete3 \
         pandas


# minimap2 binary 
RUN wget https://github.com/lh3/minimap2/releases/download/v2.21/minimap2-2.21_x64-linux.tar.bz2 && \
    tar -xf minimap2-2.21_x64-linux.tar.bz2 && \
    mv minimap2-2.21_x64-linux/minimap2 /usr/bin/

RUN apt install -y zlib1g-dev 

RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=/usr/bin/ --without-curses --disable-bz2 --disable-lzma  && \
    make && \
    make install && \
    mv samtools /usr/bin/

RUN apt install -y libjpeg-dev zlib1g-dev

RUN /usr/local/miniconda/bin/python3 -m  pip install sourmash \
    scipy \
    pysam \
    pandas \
    ete3 \
    typing \
    biopython \
    cython

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit && \
    chmod +x faSplit && \
    mv faSplit /usr/bin/

# Install conda, bowtie2, and krakenuniq
#RUN cd /usr/local/ && \
#    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
#    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
#    rm Miniconda3-latest-Linux-x86_64.sh && \
#    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
#    conda init bash && \
#    /bin/bash -c "source /root/.bashrc" && \
#    conda install -c bioconda ucsc-fasplit seqkit blast hmmer  seqtk 

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
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

ENV star_version 2.7.9a
WORKDIR /docker_main
ADD https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz /usr/bin/
RUN tar -xzf /usr/bin/${star_version}.tar.gz -C /usr/bin/
RUN cp /usr/bin/STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin

RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.95.tar.gz && \
    tar -xzf BBMap_38.95.tar.gz && \
    rm BBMap_38.95.tar.gz
RUN cp -r bbmap/* /usr/local/bin

RUN wget -q -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip/download; \ 
	unzip bowtie2.zip -d /opt/; \ 
	ln -s /opt/bowtie2-2.3.4.1-linux-x86_64/ /opt/bowtie2; \ 
	rm bowtie2.zip 
ENV PATH $PATH:/opt/bowtie2 

