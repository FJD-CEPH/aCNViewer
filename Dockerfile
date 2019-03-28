#
# Python Dockerfile
#
# https://github.com/dockerfile/python
#

# Pull base image.
FROM ubuntu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN \
  echo "tzdata tzdata/Areas select Europe\ntzdata tzdata/Zones/Europe select Berlin" > tz.txt && debconf-set-selections tz.txt

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
     x11-utils \
     && rm -rf /var/lib/apt/lists/*

# Define working directory.
WORKDIR /data

RUN \
  git clone https://github.com/FJD-CEPH/aCNViewer.git

RUN \
  python aCNViewer/code/aCNViewer.py -P installDependencies

#RUN \
#  wget http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz
RUN \
   wget https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz

RUN \
  tar xzf aCNViewer_DATA.tar.gz && rm aCNViewer_DATA.tar.gz

RUN \
  cd aCNViewer_DATA/bin/samtools-0.1.19 && make clean && make

RUN \
  cd aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_Installer && ./install -mode silent -agreeToLicense yes -destinationFolder /data/aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_ROOT

ENV LD_LIBRARY_PATH="/data/aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_ROOT/v83/runtime/glnxa64:/data/aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_ROOT/v83/bin/glnxa64:/data/aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_ROOT/v83/sys/os/glnxa64:${LD_LIBRARY_PATH}"

ENV XAPPLRESDIR=/data/aCNViewer_DATA/bin/GISTIC_2.0.23/MCR_ROOT/v83/X11/app-defaults

# Define default command.
#CMD ["python", "aCNViewer/code/aCNViewer.py"]
ENTRYPOINT ["python", "aCNViewer/code/aCNViewer.py"]
#CMD ["bash"]
