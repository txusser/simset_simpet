ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE} as base


LABEL \
  author.name="Jesús Silva"                 \
  author.email=jesus@qubiotech        \
  maintainer.email=jesus@qubiotech         \
  maintainer.url=https://scholar.google.es/citations?user=n4TyKw0AAAAJ&hl=es/  \
  source.url=https://github.com/txusser/simpet  \
  licence="MPLv2.0 (https://www.mozilla.org/en-GB/MPL/2.0/)"  \
  description="SimPET for GOLEM Decentralized Supercomputer"

RUN apt-get update \
 && apt-get install wget unzip make gcc -y

ENV LANG en_GB.UTF-8
ENV LANGUAGE en_GB:en

# This line will download the current version on simset from the original source
RUN wget -q https://github.com/txusser/simset_docker/archive/refs/heads/main.zip
RUN unzip main.zip
RUN rm main.zip

RUN sh /simset_docker-main/make_all.sh



