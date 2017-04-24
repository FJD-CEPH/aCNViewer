#
# Python Dockerfile
#
# https://github.com/dockerfile/python
#

# Pull base image.
FROM ubuntu

# Install Python.
RUN \
  apt-get update && \
  apt-get install -y \
     python \
     git \
     wget \
     python-dev \
     python-pip \
     python-virtualenv \
     r-base-core \
     r-base-dev \
     r-cran-ggplot2 \
     r-cran-gplots \
     r-cran-rcolorbrewer \
     r-cran-plotrix \
     r-bioc-limma \
     && rm -rf /var/lib/apt/lists/*

# Define working directory.
WORKDIR /data

RUN \
  git clone https://github.com/FJD-CEPH/aCNViewer.git

#RUN \
#  wget http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz
RUN \
   wget https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz

RUN \
  tar xzf aCNViewer_DATA.tar.gz

RUN \
  wget -O samtools.tar.bz2 https://sourceforge.net/projects/samtools/files/latest/download?source=files

RUN \
  cd aCNViewer_DATA/bin && tar xjf ../../samtools.tar.bz2 && make

RUN \
  python aCNViewer/code/aCNViewer.py -P installDependencies

# Define default command.
#CMD ["python", "aCNViewer/code/aCNViewer.py"]
ENTRYPOINT ["python", "aCNViewer/code/aCNViewer.py"]
#CMD ["bash"]
