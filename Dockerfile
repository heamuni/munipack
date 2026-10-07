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

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4t64 ca-certificates \
        libcfitsio10t64 libwcs8 libminpack1 libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/munipack /opt/munipack
RUN echo /opt/munipack/lib > /etc/ld.so.conf.d/munipack.conf && ldconfig
ENV PATH=/opt/munipack/bin:$PATH
WORKDIR /data
ENTRYPOINT ["munipack"]
CMD ["--help"]
