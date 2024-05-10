FROM ubuntu:23.04
ENV TAG="v1.0"


RUN apt-get update \
    && apt-get install -y procps build-essential \
    wget perl python3

###########
# ANNOVAR #
###########
ENV TOOLS=/home/TOOLS/tools
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=current
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


# INSTALL ANNOVAR
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL --wildcards *pl ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp *pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    mkdir /databases ; \
    ln -s /databases/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ;

# Download databases
# WORKDIR "$TOOLS/$TOOL_NAME/current/bin"
# RUN annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
# RUN annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad312_genome humandb/
# RUN annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20221231 humandb/
# RUN annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/

# INSTALL INTERVAR
RUN wget https://github.com/WGLab/InterVar/archive/refs/tags/v2.2.1.tar.gz && \
    tar -xvzf v2.2.1.tar.gz && \
    rm v2.2.1.tar.gz && \
    chmod a+x InterVar-2.2.1/*
RUN sed -i 's/#!\/usr\/bin\/env python/#!\/usr\/bin\/env python3/g' InterVar-2.2.1/Intervar.py
RUN mv InterVar-2.2.1/Intervar.py  $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/
RUN mv InterVar-2.2.1/intervardb $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/
WORKDIR  $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/intervardb
RUN wget https://www.omim.org/static/omim/data/mim2gene.txt
# RUN mv InterVar-2.2.1/config.ini $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/
COPY config.ini $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/

# # clean up
RUN apt-get update && apt-get upgrade -y
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# Non-root user
ARG USERNAME=docker-user
ARG USER_UID=1234
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
# Add user to root group - needed for AWS Omics(?)
RUN usermod -aG root ${USERNAME}

USER $USERNAME
ENV ANNOVAR_PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin
WORKDIR "$TOOLS/$TOOL_NAME/current/bin"

CMD ["/bin/bash"]
