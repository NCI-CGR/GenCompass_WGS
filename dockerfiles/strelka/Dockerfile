FROM ubuntu:23.04

LABEL \
  version="2.9.10" \
  description="Strelka image for use in Workflows"

# user environmental variables
ARG USERNAME=strelka-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
    
# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

# Install strelka2 dependencies
RUN apt-get update && apt-get install -y \
  bzip2 \
  g++ \
  make \
  perl-doc \
  python2 \
  rsync \
  wget \
  libncurses5-dev \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev

RUN apt-get update && apt-get upgrade -y

#HTSLIB Installation - includes bgzip and tabix
ENV HTSLIB_VESION 1.16
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VESION}/htslib-${HTSLIB_VESION}.tar.bz2 \
    && tar xvjf htslib-${HTSLIB_VESION}.tar.bz2 \
    && cd htslib-${HTSLIB_VESION} \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd / && rm -rf /tmp/htslib-${HTSLIB_VESION}

# Strelka Installation
ENV STRELKA_INSTALL_DIR /opt/strelka/
ENV STRELKA_VERSION 2.9.10

RUN wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.centos6_x86_64.tar.bz2 \
  && tar xf strelka-${STRELKA_VERSION}.centos6_x86_64.tar.bz2 \
  && rm -f strelka-${STRELKA_VERSION}.centos6_x86_64.tar.bz2 \
  && mv strelka-${STRELKA_VERSION}.centos6_x86_64 $STRELKA_INSTALL_DIR


RUN chmod a+x ${STRELKA_INSTALL_DIR}/bin/*
RUN chmod -R o+rX ${STRELKA_INSTALL_DIR}

ENV PATH="$PATH:${STRELKA_INSTALL_DIR}/bin"
# Set the default user
USER $USERNAME