FROM ubuntu:22.04
MAINTAINER timokoch@uio.no
RUN rm -f /etc/apt/apt.conf.d/docker-gzip-indexes \
  && rm -rf /var/lib/apt/lists/*

RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update \
  && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes \
  ca-certificates \
  ssh \
  vim \
  python3-dev \
  python3-pip \
  python3-venv \
  libffi-dev \
  git \
  pkg-config \
  build-essential \
  gfortran \
  cmake \
  mpi-default-dev \
  mpi-default-bin \
  libtbb-dev \
  libblas-dev \
  libsuitesparse-dev \
  libsuperlu-dev \
  libscotchparmetis-dev \
  libtrilinos-zoltan-dev \
  zlib1g-dev \
  libboost-all-dev \
  ninja-build \
  wget \
  clang \
  gmsh \
  curl \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# add CI user
RUN adduser --disabled-password --home /dumux --uid 50000 dumux

# install Python dependencies of the test suite
RUN pip3 install --upgrade pip numpy black~=22.0 pylint~=2.14.0 flake8~=6.0 matplotlib codespell~=2.2.2 fieldcompare[all]~=0.1.0

# set versions
ARG DUNE_BRANCH=releases/2.9
ARG DUMUX_BRANCH=releases/3.8

# add dune core sources
RUN mkdir -p /dune/modules
RUN cd /dune/modules && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/core/dune-common.git && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/core/dune-geometry.git && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/core/dune-grid.git && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/core/dune-localfunctions.git && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/core/dune-istl.git && \
    cd

RUN mkdir -p /dune/bin
RUN ln -s /dune/modules/dune-common/bin/dunecontrol /dune/bin/dunecontrol
ENV PATH=/dune/bin:$PATH
ENV DUNE_CONTROL_PATH=.:/dune/modules

ARG OPTS_FILE=/opts/gcc.opts
ENV DUNE_OPTS_FILE=$OPTS_FILE

# staging, extension modules and dumux
RUN cd /dune/modules && \
    git clone -b $DUNE_BRANCH --depth 1 https://gitlab.dune-project.org/extensions/dune-foamgrid.git && \
    git clone -b $DUMUX_BRANCH --depth 1 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git && \
    cd

# add the dune control option files
RUN mkdir -p /opts
COPY gcc.opts /opts/gcc.opts

# pre-build dune modules
RUN dunecontrol configure
RUN dunecontrol make all

USER dumux
WORKDIR dumux
