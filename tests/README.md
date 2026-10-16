# Tests

Munipack has two sets of tests. `make check` is the upstream test: it runs a short pipeline on one bundled real frame during the build. The tests in `tests/` exercise the installed command-line tools, one file per area (`tests/t_*.sh`).

The tests in `tests/` simulate their own frames with `munipack artificial`, so the right answer is known exactly. The tests need no network and all random numbers are seeded, so each run gives the same result. The whole set takes a few minutes.

## Running the tests

### In Docker

```bash
docker build --target test .
```

This builds the image, runs `make check`, and then runs all tests in `tests/`. A failed test fails the build. Add `--progress=plain` to see the log, and `--no-cache-filter test` to rerun the tests when nothing has changed. A plain `docker build .` skips the tests.

The final image does not contain the tests. To run them, or only some of them, in a container, tag the test stage:

```bash
docker build --target test -t munipack:test .
docker run --rm --entrypoint /opt/munipack-tests/run.sh munipack:test t_cli t_calib
```

With an image you already have, mount the test directory instead:

```bash
docker run --rm -v "$PWD/tests":/tests:ro --entrypoint /tests/run.sh <image> t_pipeline
```

### Without Docker

`tests/run.sh` tests whichever `munipack` is on your `PATH`. It needs only `bash`, `awk` and the usual coreutils.

```bash
tests/run.sh                      # all tests
tests/run.sh t_cli t_astrometry   # selected tests (.sh is optional)
KEEP=1 tests/run.sh t_gphot       # keep the working directories
```

Each test runs in its own temporary directory, which is removed after a pass and kept after a failure so you can inspect the files; the log shows its path. Each check prints one `ok:` line, a test stops at its first `FAIL:`, and `run.sh` exits with status 0 only when every test passes.

### `make check`

In a configured and built source tree (the `Dockerfile` shows the configure and build steps), run `make check` in `src/`. It finds stars, does aperture photometry, astrometry and photometric calibration on the bundled frame, and checks that the calibration constant is close to the expected value. The log is in `src/test/check.sh.log`.

## What the tests check

| Test | What it checks |
|---|---|
| `t_cli` | Help works for every action, and bad input gives a clean error (exit status 1), not a crash |
| `t_artificial` | The simulator: reproducible seeds, the detector noise model, star density, multiple processes |
| `t_calib` | Master bias, master dark and `phcorr` at low count rates |
| `t_pipeline` | The main route `find` → `aphot` → `astrometry` → `phcal`, and `psfphot` → `phcal` |
| `t_astrometry` | Astrometry from scratch on rotated frames at high declination |
| `t_astrometry_wcs` | Matching guided by the header WCS, blind matching, wrong fields, few stars, proper motions |
| `t_astrometry_roll` | Matching by a known pointing with unknown rotation, large catalogues, thread independence |
| `t_kombine` | Stacking preserves flux, positions and exposure time |
| `t_crowding` | Apertures of stars next to brighter neighbours |
| `t_phcal_filters` | Calibration in two filters gives the same result in either order |
| `t_colterm` | The colour term of one filter, and `phfotran --model mag` |
| `t_colsdss` | The colour term with real SDSS star colours |
| `t_colour2` | Two colour terms at once with real SDSS colours |
| `t_colerr` | The colour term against noisy catalogue colours: corrected for their errors, and diluted with `--color-err-correct no` |
| `t_gphot` | Growth-curve photometry, its errors and its calibration |
| `t_syserr` | The systematic error of `phcal` estimated from only a few stars |
| `t_psfempir` | Empirical PSF photometry with a PSF that varies across the frame |

## Writing a test

A test is a bash fragment `tests/t_<area>.sh`, with no shebang. `run.sh` picks up new files automatically and loads `tests/lib.sh` first, which turns on `set -e -o pipefail` and provides these helpers:

| Helper | What it does |
|---|---|
| `ok MSG`, `fail MSG` | Report a passed check, or fail the test |
| `check NAME VALUE LO HI` | Pass when VALUE is a number between LO and HI |
| `expect_fail NAME CMD...` | Pass when CMD exits with a clean error (status 1) |
| `key FILE[EXT] KEY` | The value of a header keyword |
| `nrows`, `colnum`, `column` | Number of rows, index of a column, or the values of columns of a table |
| `imstats FILE [section]` | Mean, standard deviation and number of pixels of an image |
| `imsum FILE` | A checksum of the pixels, to compare images exactly |
| `mkcat FILE N SEED [...]` | Write a synthetic catalogue of N stars with V and B magnitudes |
| `catfits`, `tabfits` | Write a catalogue or any table from text lines on stdin |
| `match` | Join two tables by sky position |
| `stats` | Median and robust standard deviation of numbers on stdin |
| `CATCOLS` | The column options that go with a `mkcat` catalogue |

Start each test with a short comment saying what it simulates. Send the output of commands that must succeed to `/dev/null`, and test every result with `check`, `expect_fail` or `fail`.
