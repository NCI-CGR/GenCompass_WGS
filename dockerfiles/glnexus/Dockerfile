FROM ubuntu:23.04

RUN apt-get update && apt-get install -y wget libjemalloc2 bcftools tabix pv
RUN apt-get update  && apt-get upgrade -y

# user environmental variables
ARG USERNAME=docker-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

# # Install GLNexus
ENV GLNEXUS_VERSION="v1.4.1"
RUN mkdir /opt/bin
RUN cd /opt/bin \
    && wget https://github.com/dnanexus-rnd/GLnexus/releases/download/${GLNEXUS_VERSION}/glnexus_cli
RUN chmod a+x /opt/bin/*
RUN chmod -R o+rX /opt
ENV PATH="$PATH:/opt/bin/"

ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
# COPY --from=builder /GLnexus/glnexus_cli /usr/local/bin/
ADD https://github.com/mlin/spVCF/releases/download/v1.3.0/spvcf /usr/local/bin/
RUN chmod +x /usr/local/bin/spvcf


ENV WORKDIR=/home/${USERNAME}
WORKDIR  ${WORKDIR}


USER $USERNAME

