# # SAMTOOLS INSTALLATION
# FROM ubuntu:23.04 AS download-samtools
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG SAMTOOLS_VERSION=1.21
# RUN curl -OL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
# RUN tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2

# FROM ubuntu:23.04 AS buildenv-samtools
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG SAMTOOLS_VERSION=1.21
# COPY --from=download-samtools /samtools-${SAMTOOLS_VERSION} /samtools-${SAMTOOLS_VERSION}
# WORKDIR /samtools-${SAMTOOLS_VERSION}
# RUN ./configure --prefix=/usr
# RUN make -j4
# RUN make install DESTDIR=/dest

# # BCFTOOLS INSTALLATION
# FROM ubuntu:23.04 AS download-bcftools
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG BCFTOOLS_VERSION=1.21
# RUN curl -OL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
# RUN tar xvjf bcftools-${BCFTOOLS_VERSION}.tar.bz2

# FROM ubuntu:23.04 AS buildenv-bcftools
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG BCFTOOLS_VERSION=1.21
# COPY --from=download-bcftools /bcftools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
# WORKDIR /bcftools-${BCFTOOLS_VERSION}
# RUN ./configure --prefix=/usr
# RUN make -j4
# RUN make install DESTDIR=/dest

# #HTSLIB Installation - includes bgzip and tabix
# FROM ubuntu:23.04 AS download-htslib
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG HTSLIB_VERSION=1.21
# RUN curl -OL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
# RUN tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2

# FROM ubuntu:23.04 AS buildenv-htslib
# RUN apt-get update \
#     && apt-get install -y procps build-essential \
#     wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
# ARG HTSLIB_VERSION=1.21
# COPY --from=download-htslib /htslib-${HTSLIB_VERSION} /htslib-${HTSLIB_VERSION}
# WORKDIR /htslib-${HTSLIB_VERSION}
# RUN ./configure --prefix=/usr
# RUN make -j4
# RUN make install DESTDIR=/dest




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
    bcftools==1.21



RUN conda install -c conda-forge \
    numpy==2.0.0 \
    pandas==2.2.2 \
    xlsxwriter==3.1.9 \
    openpyxl==3.1.4 \
    python-dateutil==2.9.0 \
    cryptography==42.0.8 \
    click==8.1.7 \
    && conda clean -afy


# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY scripts/* /opt/bin/
RUN chmod a+x /opt/bin/*
RUN chmod -R o+rX /opt


# RUN  apt remove -y python3-pip
RUN pip uninstall -y pip
# ENV WORKDIR=/home/${USERNAME}
# WORKDIR  ${WORKDIR}

ENV PATH="$PATH:/opt/bin/"

USER $USERNAME

