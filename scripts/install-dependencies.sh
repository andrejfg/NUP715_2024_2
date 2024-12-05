#!/bin/bash
apt-get update 
apt-get -y install default-jre autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev libperl-dev libgsl0-dev tabix

chmod +x tools/prinseq-lite
chmod +x tools/fastqc/fastqc