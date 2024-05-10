# SAMTOOLS INSTALLATION
FROM ubuntu:23.04 AS download-samtools
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG SAMTOOLS_VERSION=1.16.1
RUN curl -OL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2

FROM ubuntu:23.04 AS buildenv-samtools
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG SAMTOOLS_VERSION=1.16.1
COPY --from=download-samtools /samtools-${SAMTOOLS_VERSION} /samtools-${SAMTOOLS_VERSION}
WORKDIR /samtools-${SAMTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

# BCFTOOLS INSTALLATION
FROM ubuntu:23.04 AS download-bcftools
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG BCFTOOLS_VERSION=1.16
RUN curl -OL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar xvjf bcftools-${BCFTOOLS_VERSION}.tar.bz2

FROM ubuntu:23.04 AS buildenv-bcftools
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG BCFTOOLS_VERSION=1.16
COPY --from=download-bcftools /bcftools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
WORKDIR /bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

#HTSLIB Installation - includes bgzip and tabix
FROM ubuntu:23.04 AS download-htslib
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG HTSLIB_VERSION=1.16
RUN curl -OL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
RUN tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2

FROM ubuntu:23.04 AS buildenv-htslib
RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl
ARG HTSLIB_VERSION=1.16
COPY --from=download-htslib /htslib-${HTSLIB_VERSION} /htslib-${HTSLIB_VERSION}
WORKDIR /htslib-${HTSLIB_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest


FROM ubuntu:23.04

RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget zlib1g bzip2 make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev parallel curl

COPY --from=buildenv-bcftools /dest /
COPY --from=buildenv-htslib /dest /
COPY --from=buildenv-samtools /dest /

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
    snpeff==5.1 \
     snpsift==5.1


RUN conda install -c conda-forge openjdk==20.0.0
RUN snpEff download GRCh38.105

# RUN  apt remove -y python3-pip
RUN pip uninstall -y pip
# ENV WORKDIR=/home/${USERNAME}
# WORKDIR  ${WORKDIR}

ENV PATH="$PATH:/opt/bin/"

USER $USERNAME

