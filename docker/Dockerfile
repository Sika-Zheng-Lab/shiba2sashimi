##################################################################
# Dockerfile to build container images for shiba2sashimi (v0.1.0)
# Based on python 3.9-bullseye
##################################################################

FROM python:3.9-bullseye

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

# Install Arial font (msttcorefonts)
RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
    fontconfig && \
    wget -q http://ftp.de.debian.org/debian/pool/contrib/m/msttcorefonts/ttf-mscorefonts-installer_3.8_all.deb && \
    apt-get -qq -y install --no-install-recommends \
    cabextract && \
    dpkg -i ttf-mscorefonts-installer_3.8_all.deb || true && \
    apt-get -qq -y --fix-broken install && \
    fc-cache -fv && \
    rm -rf ttf-mscorefonts-installer_3.8_all.deb

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
