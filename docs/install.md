# Installation

The recommended way to run Munipack is the Docker image built from this
repository. It contains the command-line tools only (no GUI) on Ubuntu 24.04.
To install natively instead, see {ref}`building-from-source`.

## Docker image

### Build

From the top of the repository:

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
```

The first build takes about 10 minutes. Later builds reuse the cached
dependency layers. The image is about 450 MB, mostly Open MPI and the
packages it pulls in.

The Dockerfile has these stages:

| Stage | Contents |
|---|---|
| `build` | Ubuntu 24.04 with the compilers and `-dev` packages. Builds wxWidgets 3.3.3 (non-GUI part only, checksum verified), OakLeaf from `deps/oakleaf`, then Munipack from `src/`, installed in `/opt/munipack`. Runs the upstream test suite (`make check`); if it fails, the build stops and prints the test log. |
| `runtime` | A clean Ubuntu 24.04 with only the runtime libraries and `/opt/munipack` copied from `build`. Entry point `munipack`, working directory `/data`. |
| `test` | `runtime` plus `tests/`. Runs `tests/run.sh` against the installed tools. |
| *(final)* | The same as `runtime`. A plain `docker build` produces this one. |

A plain build skips the `test` stage, because BuildKit only builds the stages
the target needs. To run the tests as well:

```bash
docker build --target test .
```

The tests use frames simulated by `artificial`, so they need no network. They
take 2–3 minutes. {doc}`tests` describes each one.

```{note}
wxWidgets is downloaded from GitHub during the build, so the build needs
network access. Upstream revision hg1653 uses the wxWidgets 3.3 API, but
Ubuntu 24.04 ships 3.2, which is why the build compiles it from source.
```

### Running

Mount the directory with your data at `/data` and run as your own user, so
the output files belong to you and not to root:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes --help
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes find frame.fits
```

The entry point is `munipack`. Its first argument is the action (`find`,
`aphot`, `astrometry`, ...) and `munipack <action> --help` lists that
action's options. Paths are relative to `/data`, so files outside the mounted
directory are not visible. Mount another directory with a second `-v` if you
need one.

A shell alias saves typing:

```bash
alias munipack='docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes'
munipack aphot frame.fits
```

`cone` and `sesame` query Virtual Observatory services over the network. This
works with Docker's default network; a `--network none` container cannot
reach them.

## Parallel processing

Munipack has two independent kinds of parallelism.

### Processes: `--mproc`

`artificial` and `astrometry` can spread their frames over `--mproc N`
processes. These are Fortran coarrays, run on Open MPI through OpenCoarrays.

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
  artificial -c cone.fits --count 8 --seed 3 --mproc 4
```

- The default is one process, set in the image by `MUNIPACK_CAF_MPROC=1`.
  With one process, Munipack runs the programs built for a single process,
  without MPI. With `--mproc N` greater than 1, it runs their parallel
  variants (`artificial-caf`, `astrometry-caf`) under `mpirun -n N`.
- N is capped at the number of logical CPUs.
- The results do not depend on N. Frames made with `artificial --seed` are
  bit-identical for any number of processes, and `astrometry` gives the same
  solutions.
- Each process handles whole frames, so `--mproc` only helps when there are
  several frames.
- `--Xproc "<options>"` (or the `MUNIPACK_CAF_XPROC` variable) passes extra
  options to `mpirun`. `MUNIPACK_CAF_RUN` changes the launcher.
- The image sets the Open MPI variables a container needs (see the
  `ENV` block of the Dockerfile). They turn off UCX one-sided communication,
  allow `mpirun` as root, allow more processes than cores, and suppress help
  banners.

### Threads: OpenMP

The blind match of `astrometry` runs in OpenMP threads within each process.
It is the search used without a usable WCS or pointing, after a failed
guided match, or with `--blind`.

- `OMP_NUM_THREADS` sets the number of threads. When it is unset, Munipack
  uses the CPUs available to the process divided by `--mproc`. The available
  CPUs take account of the affinity mask and the cgroup CPU quota, so
  `docker run --cpus 2` gives two threads.
- The solution is the same for any number of threads. `OMP_NUM_THREADS=1`
  reproduces the single-threaded behaviour.
- The threads share out the pairs of chunks (of the frame and of the
  catalogue). One pair is never split between threads, so a single slow
  pair sets a floor on the run time.

```bash
docker run --rm --cpus 4 -u $(id -u):$(id -g) -v "$PWD":/data \
  -e OMP_NUM_THREADS=4 munipack:0.6.2-hg1653-fixes astrometry -c cat.fits frame.fits
```

(building-from-source)=

## Building from source

The Dockerfile is the reference recipe. In outline, on Ubuntu 24.04:

1. Install the build dependencies:

   ```bash
   sudo apt-get install build-essential gfortran pkg-config autoconf automake libtool \
       ca-certificates curl bzip2 libcurl4-openssl-dev \
       libcfitsio-dev wcslib-dev minpack-dev \
       libcoarrays-openmpi-dev libopenmpi-dev openmpi-bin
   ```

   The coarray packages are only needed for `--mproc`. Use Open MPI: with
   Ubuntu's MPICH 4.2.0, every process runs alone as image 1 of 1. Do not
   install `libcoarrays-mpich-dev` as well, because `configure` takes the
   MPICH coarray library first when it finds one.

2. Set the environment, so that each of the next steps finds the previous
   ones:

   ```bash
   export PREFIX=/opt/munipack
   export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig LDFLAGS=-L$PREFIX/lib \
          LD_LIBRARY_PATH=$PREFIX/lib PATH=$PREFIX/bin:$PATH \
          ACLOCAL_PATH=$PREFIX/share/aclocal
   ```

3. Build wxWidgets 3.3 (3.3.3 in the image) without the GUI, in its
   unpacked source:

   ```bash
   curl -fsSLO https://github.com/wxWidgets/wxWidgets/releases/download/v3.3.3/wxWidgets-3.3.3.tar.bz2
   tar xjf wxWidgets-3.3.3.tar.bz2 && cd wxWidgets-3.3.3
   ./configure --prefix=$PREFIX --disable-gui \
       --with-zlib=builtin --with-regex=builtin --with-expat=builtin \
       --enable-webrequest --without-libmspack
   make -j"$(nproc)" && make install
   ```

4. Build OakLeaf in `deps/oakleaf`. Ubuntu's `liboakleaf` is too old and
   has no pkg-config file.

   ```bash
   autoreconf -i && ./configure --prefix=$PREFIX && make -j"$(nproc)" && make install
   ```

5. Build Munipack in `src/`. The build has to be in-tree, because the
   documentation Makefile does not support VPATH builds.

   ```bash
   ./bootstrap --dev 0.6.2+hg1653
   autoreconf -i
   ./configure --prefix=$PREFIX --disable-gui FCFLAGS="-O2 -I$PREFIX/include"
   make -j"$(nproc)" && make install
   make check
   ```

At the end, `configure` reports whether it found coarray support and the
OpenMP flags. Two options change these:

- `--disable-caf` builds without coarrays, so only one process is available.
- `--disable-openmp` builds the blind match with a single thread.

```{note}
In a native build, `--mproc` defaults to 3 processes when the MPI coarray
variant was built. Set `MUNIPACK_CAF_MPROC=1` to make one process the
default, as the image does. Outside a container, you normally do not need
the `OMPI_ALLOW_RUN_AS_ROOT*` variables. The image also sets
`OMPI_MCA_osc=^ucx`, because UCX one-sided communication failed there in
`MPI_Win_create`; set it too if `mpirun` fails the same way.
```

To run this fork's tests against the `munipack` on `PATH`:

```bash
tests/run.sh                 # all tests
tests/run.sh t_cli t_calib   # selected tests
KEEP=1 tests/run.sh t_pipeline   # keep the working directories
```
