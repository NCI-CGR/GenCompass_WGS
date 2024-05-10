FROM ubuntu:23.04

# Non-root user
ARG USERNAME=docker-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

# install stuff with apt-get
RUN apt-get update && apt-get install -y procps build-essential  wget \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get upgrade -y

WORKDIR /usr/src

#download fastp binary
ENV FASTP_VERSION="0.23.4"
RUN mkdir fastp \
    && cd fastp \
    && wget http://opengene.org/fastp/fastp.${FASTP_VERSION} \
    && mv fastp.${FASTP_VERSION} fastp \
    && chmod a+x ./fastp

ENV PATH="/usr/src/fastp/:$PATH"

# RUN  apt remove -y python3-pip
# copy bash/python scripts specific to your pipeline
# COPY scripts/* /usr/local/bin/

# USER $USERNAME



