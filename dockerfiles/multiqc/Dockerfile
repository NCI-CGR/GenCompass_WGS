FROM ubuntu:23.04

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
# parallel to run GNU parallel processes
# pip to install python packages
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl

LABEL description="GenCompass custom scripts"
ARG ENV_NAME="gencompass-base"


# user environmental variables
ARG USERNAME=gencompass-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

RUN apt-get update && apt-get upgrade -y
RUN apt-get clean -y && rm -rf /var/lib/apt/lists/*
# install dependency required by bcftools
# RUN apt-get install --yes wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev

# INSTALL CONDA 
ENV CONDAPATH=/usr/local/bin/miniconda3
ENV PATH=$CONDAPATH/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh \
    && bash ./Miniconda3*.sh -b -p $CONDAPATH \
    && rm Miniconda3*.sh



RUN conda install -c bioconda \
    multiqc==1.15 \
    && conda clean -afy


# RUN  apt remove -y python3-pip
RUN pip uninstall -y pip
# ENV WORKDIR=/home/${USERNAME}
# WORKDIR  ${WORKDIR}

# USER $USERNAME

