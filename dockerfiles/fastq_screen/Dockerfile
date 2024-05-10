FROM ubuntu:23.04

RUN apt-get update

RUN apt-get install -y build-essential zlib1g zlib1g-dev make gcc g++ wget \
   bowtie2 bwa libgd-perl pkg-config libgd-dev cmake cpanminus

RUN apt-get update && apt-get upgrade -y

RUN cpanm GD::Graph::bars

ENV FASTQ_SCREEN_VERSION="0.15.3"
ADD https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v${FASTQ_SCREEN_VERSION}.tar.gz /opt
WORKDIR /opt
RUN tar xzvf v${FASTQ_SCREEN_VERSION}.tar.gz
WORKDIR /opt/FastQ-Screen-${FASTQ_SCREEN_VERSION}
RUN chmod 755 fastq_screen
ENV PATH="/opt/FastQ-Screen-${FASTQ_SCREEN_VERSION}:$PATH"

RUN apt-get clean

# # user environmental variables
# ARG USERNAME=gencompass-user
# ARG USER_UID=1122
# ARG USER_GID=$USER_UID

# # Create the user
# RUN groupadd --gid $USER_GID $USERNAME \
#     && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# # Add user to root group - needed for AWS Omics(?)
# RUN usermod -aG root ${USERNAME}

# USER $USERNAME
