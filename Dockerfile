# Build Munipack command-line tools (no GUI) from ./src: an upstream
# Mercurial snapshot plus local changes.
FROM ubuntu:24.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran pkg-config autoconf automake libtool \
        ca-certificates curl bzip2 \
        libcurl4-openssl-dev \
        libcfitsio-dev wcslib-dev minpack-dev \
    && rm -rf /var/lib/apt/lists/*

ENV PREFIX=/opt/munipack
ENV PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig \
    LDFLAGS=-L$PREFIX/lib \
    LD_LIBRARY_PATH=$PREFIX/lib \
    PATH=$PREFIX/bin:$PATH \
    ACLOCAL_PATH=$PREFIX/share/aclocal

# The snapshot uses wxWidgets 3.3 API (Ubuntu ships 3.2). Only the non-GUI
# libraries (base, net, xml) are built: the CLI and VO tools need nothing else.
ARG WX_VERSION=3.3.3
ARG WX_SHA256=81b09d6dd9f1ed9301f8c55a968a488d0491f264dc2bab19a7e407ac67009482
WORKDIR /build/wx
RUN curl -fsSLO https://github.com/wxWidgets/wxWidgets/releases/download/v${WX_VERSION}/wxWidgets-${WX_VERSION}.tar.bz2 \
    && echo "${WX_SHA256}  wxWidgets-${WX_VERSION}.tar.bz2" | sha256sum -c - \
    && tar xjf wxWidgets-${WX_VERSION}.tar.bz2 --strip-components=1 \
    && ./configure --prefix=$PREFIX --disable-gui \
        --with-zlib=builtin --with-regex=builtin --with-expat=builtin \
        --enable-webrequest --without-libmspack \
    && make -j"$(nproc)" \
    && make install

# Ubuntu's liboakleaf is too old, and the Munipack snapshot needs OakLeaf's
# pkg-config file, which only the upstream hg snapshot in ./deps/oakleaf has.
COPY deps/oakleaf /build/oakleaf
WORKDIR /build/oakleaf
RUN autoreconf -i \
    && ./configure --prefix=$PREFIX \
    && make -j"$(nproc)" \
    && make install

# Coarray Fortran on Open MPI (OpenCoarrays) for the parallel processes
# (--mproc) of artificial and astrometry; configure finds it by pkg-config.
# Ubuntu's MPICH 4.2.0 does not work: under its mpirun every process runs
# alone as image 1 of 1.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcoarrays-openmpi-dev libopenmpi-dev openmpi-bin \
    && rm -rf /var/lib/apt/lists/*

# Version string the upstream bootstrap would derive from `hg tags`.
ARG MUNIPACK_VERSION=0.6.2+hg1653
COPY src /build/src
WORKDIR /build/src
# Built in-tree: doc/Makefile.am does not support out-of-tree (VPATH) builds.
RUN ./bootstrap --dev "${MUNIPACK_VERSION}" \
    && autoreconf -i \
    && ./configure --prefix=$PREFIX --disable-gui FCFLAGS="-O2 -I$PREFIX/include" \
    && make -j"$(nproc)" \
    && make install
# Upstream test suite (test/check.sh: find, aphot, astrometry, phcal on a
# bundled frame); a failure stops the build and prints the test log. It runs
# with the Open MPI settings of the final image below.
RUN MUNIPACK_CAF_MPROC=1 OMPI_MCA_osc=^ucx \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 OMPI_MCA_orte_execute_quiet=1 \
    make check || { cat test/check.sh.log; exit 1; }

FROM ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive
# Open MPI depends on UCX, which depends on the AMD GPU runtime (with LLVM,
# 150 MB), and on the hwloc plugins (libxml2 with ICU, X11; 45 MB). Neither
# is used (processes on one host, no GPU: UCX and hwloc load them as
# optional plugins), so empty placeholder packages stand in for them.
RUN for p in libamdhip64-5 libhwloc-plugins; do \
        mkdir -p /tmp/placeholder/$p/DEBIAN \
        && printf 'Package: %s\nVersion: 99\nArchitecture: all\nMaintainer: Munipack Docker build\nDescription: empty placeholder, not used by Munipack\n' $p \
           > /tmp/placeholder/$p/DEBIAN/control \
        && dpkg-deb --build /tmp/placeholder/$p /tmp/placeholder/$p.deb \
        && dpkg -i /tmp/placeholder/$p.deb || exit 1; \
    done \
    && rm -rf /tmp/placeholder
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4t64 ca-certificates \
        libcfitsio10t64 libwcs8 libminpack1 libgfortran5 \
        libcaf-openmpi-3t64 openmpi-bin \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/munipack /opt/munipack
RUN echo /opt/munipack/lib > /etc/ld.so.conf.d/munipack.conf && ldconfig
ENV PATH=/opt/munipack/bin:$PATH
# One process unless --mproc N is given (MUNIPACK_CAF_MPROC is the default):
# it runs the programs built for one process, as before, without MPI;
# --mproc N > 1 runs artificial-caf, astrometry-caf by mpirun. Open MPI:
# one-sided communication without UCX (UCX failed in MPI_Win_create with one
# process), mpirun as root (the default user of the container), more
# processes than cores (--mproc is limited by the logical CPUs, Open MPI
# counts the cores), no help banners (MPI_ABORT on STOP with a message).
ENV MUNIPACK_CAF_MPROC=1 \
    OMPI_MCA_osc=^ucx \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 \
    OMPI_MCA_orte_execute_quiet=1
WORKDIR /data
ENTRYPOINT ["munipack"]
CMD ["--help"]
