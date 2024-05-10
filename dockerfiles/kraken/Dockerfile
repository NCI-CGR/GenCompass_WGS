# ccbr_kraken_v2.1.3 Dockerfile (v0.0.1)
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



RUN mkdir -p /data2
RUN mkdir -p /opt2

# Update apt-get before downloading packages 
RUN apt-get update && \
    apt-get upgrade -y 

# Install dependencies for latest versions of Kraken (v2.1.3) and Krona (v.2.8.1)
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
	build-essential \
	wget \
	curl \
	zlib1g-dev \
	make \
	python3 \
	perl \
	g++ \
	locales

RUN apt-get update && apt-get upgrade -y
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt2

# Get latest kraken release (v2.1.3): https://github.com/DerrickWood/kraken2/releases 
ENV KRAKEN_VERSION="2.1.3"
RUN wget https://github.com/DerrickWood/kraken2/archive/v${KRAKEN_VERSION}.tar.gz && \
	tar -xzf v${KRAKEN_VERSION}.tar.gz && \
	rm v${KRAKEN_VERSION}.tar.gz && \
    mv kraken2-${KRAKEN_VERSION} kraken2 && \
	cd kraken2 && \
 	./install_kraken2.sh . && \
	chmod a+rX /opt2/kraken2/*

# Add kraken2 to PATH
ENV PATH="$PATH:/opt2/kraken2"

WORKDIR /opt2

# Get latest krona release (v2.7.1): https://github.com/marbl/Krona/releases
# Build Krona Taxonomy Databases
# ktImportTaxonomy relies on NCBI taxonomy: https://github.com/marbl/Krona/wiki/Installing 
ENV KRONA_VERSION="2.8.1"
RUN wget https://github.com/marbl/Krona/archive/v${KRONA_VERSION}.tar.gz && \
	tar -xzf v${KRONA_VERSION}.tar.gz && \
	rm v${KRONA_VERSION}.tar.gz && \ 
	cd Krona-${KRONA_VERSION}/KronaTools && \
	./install.pl  && \
	mkdir -p /opt2/Krona-${KRONA_VERSION}/KronaTools/taxonomy && \
	./updateTaxonomy.sh

# Set environment variable(s)
# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Dockerfile to build image located at /opt2/Dockerfile
COPY Dockerfile /opt2/
WORKDIR /data2

# USER $USERNAME

