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
         typing


# # conda stuff 
# ENV PATH="/root/miniconda3/bin:${PATH}"
# ARG PATH="/root/miniconda3/bin:${PATH}"

# #RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# RUN wget \
#     https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#     && mkdir /root/.conda \
#     && bash Miniconda3-latest-Linux-x86_64.sh -b \
#     && rm -f Miniconda3-latest-Linux-x86_64.sh 

# # updating R 
# #RUN apt install -y --no-install-recommends software-properties-common dirmngr && \
# #    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
# #    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
# #    apt -y install r-base



# RUN apt-get -y update && \
#     apt-get install -y default-jdk
# RUN conda install -c bioconda metalign

