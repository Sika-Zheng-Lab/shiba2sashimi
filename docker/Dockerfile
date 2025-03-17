##################################################################
# Dockerfile to build container images for shiba2sashimi (v0.1.0)
# Based on python 3.7-buster
##################################################################

FROM python:3.7-buster

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
    build-essential gcc g++ make cmake libtool texinfo dpkg-dev pkg-config \
    libgsl-dev wget locate less vim zlib1g-dev bzip2 lzma curl \
    libboost-dev libcurl4-openssl-dev libboost-all-dev libbz2-dev liblzma-dev \
    libpcre3 libpcre3-dev

# Install shiba2sashimi
RUN pip install shiba2sashimi==0.1.0

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
