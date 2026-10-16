# Installation

The easiest way to run Munipack is the Docker image built from this
repository. It contains the command-line tools (no GUI) on Ubuntu 24.04, so
you do not have to build any dependencies yourself. If you prefer a native
install, see {ref}`building-from-source`.

## Docker image

### Build

From the top of the repository:

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
```

The first build takes a while (around 10 minutes) and needs network access,
because it downloads wxWidgets. Later builds are much faster. The build also
runs the upstream test suite and stops if it fails.

To run this fork's own tests as well, build the `test` target. They use
simulated frames and need no network; {doc}`tests` describes them.

```bash
docker build --target test .
```

### Run

Mount your data directory at `/data` and run as your own user, so the output
files belong to you rather than to root:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes find frame.fits
```

The first argument is the action (`find`, `aphot`, `astrometry`, ...), and
`munipack <action> --help` lists its options. Paths are relative to `/data`,
so only files in the mounted directory are visible; add another `-v` if you
need more.

An alias makes this feel like a normal command:

```bash
alias munipack='docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes'
munipack aphot frame.fits
```

`cone` and `sesame` query online catalogues, so they need the container's
network (the Docker default works; `--network none` does not).

## Running in parallel

`artificial` and `astrometry` can process several frames at once with
`--mproc N`. Each process takes whole frames, so this only helps when you
have more than one frame. The results are the same for any N, and N is
capped at the number of CPUs. `--Xproc "<options>"` passes extra options to
`mpirun`.

```bash
munipack artificial -c cone.fits --count 8 --seed 3 --mproc 4
```

The blind match in `astrometry` also uses threads. Set the number with
`OMP_NUM_THREADS`. By default Munipack uses the CPUs the container is allowed
(`docker run --cpus 2` gives two threads), shared between the `--mproc`
processes. The solution does not depend on the number of threads.

```bash
docker run --rm --cpus 4 -u $(id -u):$(id -g) -v "$PWD":/data \
  -e OMP_NUM_THREADS=4 munipack:0.6.2-hg1653-fixes astrometry -c cat.fits frame.fits
```

(building-from-source)=

## Building from source

The `Dockerfile` is the reference recipe. On Ubuntu 24.04 the steps are:

1. Install the build dependencies:

   ```bash
   sudo apt-get install build-essential gfortran pkg-config autoconf automake libtool \
       ca-certificates curl bzip2 libcurl4-openssl-dev \
       libcfitsio-dev wcslib-dev minpack-dev \
       libcoarrays-openmpi-dev libopenmpi-dev openmpi-bin
   ```

   The coarray and MPI packages are only needed for `--mproc`. Use the Open
   MPI ones: parallel runs do not work with Ubuntu's MPICH, and `configure`
   picks `libcoarrays-mpich-dev` first if it is installed too.

2. Set the environment so that each step finds the previous ones:

   ```bash
   export PREFIX=/opt/munipack
   export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig LDFLAGS=-L$PREFIX/lib \
          LD_LIBRARY_PATH=$PREFIX/lib PATH=$PREFIX/bin:$PATH \
          ACLOCAL_PATH=$PREFIX/share/aclocal
   ```

3. Build wxWidgets 3.3 without the GUI (Ubuntu ships only 3.2):

   ```bash
   curl -fsSLO https://github.com/wxWidgets/wxWidgets/releases/download/v3.3.3/wxWidgets-3.3.3.tar.bz2
   tar xjf wxWidgets-3.3.3.tar.bz2 && cd wxWidgets-3.3.3
   ./configure --prefix=$PREFIX --disable-gui \
       --with-zlib=builtin --with-regex=builtin --with-expat=builtin \
       --enable-webrequest --without-libmspack
   make -j"$(nproc)" && make install
   ```

4. Build OakLeaf in `deps/oakleaf` (Ubuntu's package is too old):

   ```bash
   autoreconf -i && ./configure --prefix=$PREFIX && make -j"$(nproc)" && make install
   ```

5. Build Munipack in `src/`, inside the source tree:

   ```bash
   ./bootstrap --dev 0.6.2+hg1653
   autoreconf -i
   ./configure --prefix=$PREFIX --disable-gui FCFLAGS="-O2 -I$PREFIX/include"
   make -j"$(nproc)" && make install
   make check
   ```

   Add `--disable-caf` to build without parallel processes, or
   `--disable-openmp` to build without threads.

```{note}
A native build with parallel support uses 3 processes by default. Set
`MUNIPACK_CAF_MPROC=1` to make it one, as in the Docker image. If `mpirun`
fails in `MPI_Win_create`, also set `OMPI_MCA_osc=^ucx`.
```

To run this fork's tests against the `munipack` on your `PATH`:

```bash
tests/run.sh                    # all tests
tests/run.sh t_cli t_calib      # selected tests
KEEP=1 tests/run.sh t_pipeline  # keep the working directories
```
