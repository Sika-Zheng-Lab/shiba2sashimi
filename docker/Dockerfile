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
    wget less vim curl

# Upgrade pip
RUN pip install --upgrade pip

# Install shiba2sashimi
RUN git clone https://github.com/Sika-Zheng-Lab/shiba2sashimi.git -b v0.1.0 && \
    cd shiba2sashimi && \
    pip install .

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
