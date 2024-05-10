#!/bin/bash

set -xeu -o pipefail

docker build -t informationsea/snpeff:5.1-GRCh38.p14 .
docker build --build-arg DB_VERSION=GRCh37.p13 -t informationsea/snpeff:5.1-GRCh37.p13 .

for one in informationsea/snpeff:5.1-GRCh38.p14 informationsea/snpeff:5.1-GRCh37.p13; do
    docker run -v /var/run/docker.sock:/var/run/docker.sock \
        -v $PWD/out:/output \
        --privileged -t --rm \
        quay.io/singularity/docker2singularity \
        --name $(echo ${one}|sed -e 's/[^a-zA-Z0-9.\-]/_/g') \
        ${one}
done