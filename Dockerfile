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
  apt-get install -y python python-dev python-pip python-virtualenv git r-base-core r-base-dev r-cran-ggplot2 r-cran-plotrix & \
  rm -rf /var/lib/apt/lists/*

# Define working directory.
WORKDIR /data

RUN \
  git clone https://github.com/FJD-CEPH/aCNViewer.git



# Define default command.
CMD ["python", "aCNViewer/code/aCNViewer.py"]
