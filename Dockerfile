FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \ 
    cmake \
    build-essential \
#    libeigen3-dev \
    gdb

ENTRYPOINT ["bash"]
