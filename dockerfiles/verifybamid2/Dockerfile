# In this image we do not include the R-based plot scripts.
FROM ubuntu:23.04 as build

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive
ENV LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

# apt update and install build requirements
RUN apt-get update \
     && apt-get install -y \
        g++ \
        cmake \
        git \
        wget \
        libbz2-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        liblzma-dev \
        bzip2

RUN apt-get update && apt-get upgrade -y


# Compile htslib
# WORKDIR /deps
# RUN wget -q https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
#     && tar -xf htslib-1.11.tar.bz2 \
#     && mv htslib-1.11 htslib
# WORKDIR /deps/htslib
# RUN autoheader; autoconf; ./configure --prefix=/usr/local/ \
#     && make && make install

# #HTSLIB Installation - includes bgzip and tabix
ENV HTSLIB_VESION 1.16
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VESION}/htslib-${HTSLIB_VESION}.tar.bz2 \
    && tar xvjf htslib-${HTSLIB_VESION}.tar.bz2 \
    && cd htslib-${HTSLIB_VESION} \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd / && rm -rf /tmp/htslib-${HTSLIB_VESION}


# Compile VerifyBamID. Version 2.0.1
WORKDIR /
# RUN git clone --depth 1 --branch 2.0.1 git://github.com/Griffan/VerifyBamID.git
# WORKDIR /VerifyBamID/build
# RUN  cmake .. \
#      && make \
#      && make test
RUN wget https://github.com/Griffan/VerifyBamID/releases/download/2.0.1/VerifyBamID.Linux.x86-64.gz \
    && gzip -d VerifyBamID.Linux.x86-64.gz
RUN mv VerifyBamID.Linux.x86-64 verifybamid2
RUN chmod a+x verifybamid2





# Final image
FROM ubuntu:23.04
# ENV DEBIAN_FRONTEND noninteractive

RUN apt update && apt-get upgrade && apt install -y \
        libgomp1 \
        libcurl4-openssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /verifybamid2 /usr/local/bin/
ENV PATH="$PATH:/usr/local/bin/"
# COPY --from=build /usr/local/lib/ /usr/local/lib/

# Non-root user
ARG USERNAME=docker-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

# CMD ["VerifyBamID"]


# USER $USERNAME

