FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages
# TODO(akhil): May need to set this up to use specific packages instead of the
# apt-repo
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    build-essential \
    gdb \
    libboost-all-dev \
    libeigen3-dev \
    libgtest-dev \
    ninja

RUN mkdir -p /usr/filter
ENV CMAKE_GENERATOR=Ninja

COPY . /usr/filter

ENTRYPOINT ["bash"]
