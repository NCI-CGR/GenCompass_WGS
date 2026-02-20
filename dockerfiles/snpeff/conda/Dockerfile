FROM ubuntu:23.10

RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl unzip zip

# COPY --from=buildenv-bcftools /dest /
# COPY --from=buildenv-htslib /dest /
# COPY --from=buildenv-samtools /dest /

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
# Install mamba for faster installation in the subsequent step
# Install r-base for being able to run the install.R script

# RUN conda install -c conda-forge mamba r-base -y

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
# parallel to run GNU parallel processes
# pip to install python packages

RUN apt-get update && apt-get upgrade -y
RUN apt-get clean -y && rm -rf /var/lib/apt/lists/*

# INSTALL CONDA 
ENV CONDAPATH=/usr/local/bin/miniconda3
ENV PATH=$CONDAPATH/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py312_24.5.0-0-Linux-x86_64.sh \
    && bash ./Miniconda3*.sh -b -p $CONDAPATH \
    && rm Miniconda3*.sh


# INSTALL AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install


RUN conda install -c bioconda -c conda-forge \
    bedtools==2.31.1 \
    pysam==0.22.1 \
    htslib==1.21 \
    samtools==1.21 \
    bcftools==1.21 \
    snpeff==5.2


RUN conda install -c conda-forge \
    numpy==2.0.0 \
    pandas==2.2.2 \
    openpyxl==3.1.4 \
    python-dateutil==2.9.0 \
    cryptography==42.0.8 \
    && conda clean -afy

# Copy additional scripts from bin and add to PATH
# RUN mkdir /opt/bin
# COPY scripts/* /opt/bin/
# RUN chmod a+x /opt/bin/*
# RUN chmod -R o+rX /opt


# RUN  apt remove -y python3-pip
RUN pip uninstall -y pip
# ENV WORKDIR=/home/${USERNAME}
# WORKDIR  ${WORKDIR}

ENV PATH="$PATH:/opt/bin/"

USER $USERNAME

