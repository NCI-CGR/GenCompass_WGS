FROM --platform=linux/amd64 ubuntu:23.04


# user environmental variables
ARG USERNAME=gencompass-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}


RUN apt-get update \
    && apt-get install -y procps build-essential wget 

RUN apt-get update && apt-get upgrade -y

RUN mkdir /opt/bin

ENV SOMALIER_VESION v0.2.16
RUN cd /opt/bin \
    && wget https://github.com/brentp/somalier/releases/download/${SOMALIER_VESION}/somalier

RUN chgrp -R +${USER_GID} /opt
RUN chmod a+x /opt/bin/*
RUN chmod -R o+rX /opt


ENV WORKDIR=/home/${USERNAME}
WORKDIR  ${WORKDIR}

ENV PATH="$PATH:/opt/bin/"
# USER $USERNAME
