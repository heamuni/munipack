# Tests

Munipack has two test suites:

- **`make check`**: the upstream check in `src/test/check.sh`. It runs on one bundled real frame during the build.
- **`tests/`**: tests of the installed command-line tools, one file per area (`tests/t_*.sh`).

Each test in `tests/` simulates its frames with `munipack artificial`, mostly from a synthetic catalogue and in two cases from real SDSS stars embedded in the test. The truth is therefore known exactly. The tests need no network, and every random number is seeded, so a run is reproducible. The whole set takes one to three minutes, depending on the machine.

## Running the tests

### In Docker

```bash
docker build --target test .
```

This builds the image, runs `make check` inside the build stage, then runs `tests/run.sh` against the installed tools in the `test` stage. A failed test fails the build. A plain `docker build .` skips the `test` stage.

- **Seeing the log:** add `--progress=plain`.
- **Running the tests again unchanged:** the stage is cached. Add `--no-cache-filter test` to force it to run again.

The final image does not contain the tests. To run them, or only some of them, in a container, tag the `test` stage and override the entrypoint:

```bash
docker build --target test -t munipack:test .
docker run --rm --entrypoint /opt/munipack-tests/run.sh munipack:test              # all tests
docker run --rm --entrypoint /opt/munipack-tests/run.sh munipack:test t_cli t_calib
```

The same works with the final image if you mount the test directory:

```bash
docker run --rm -v "$PWD/tests":/tests:ro --entrypoint /tests/run.sh \
    munipack:0.6.2-hg1653-fixes t_pipeline
```

### Without Docker

`tests/run.sh` tests the `munipack` found on `PATH`. It needs `bash`, `awk` (gawk or mawk) and the usual coreutils (`md5sum`, `mktemp`, `sort`, `paste`, `seq`).

```bash
tests/run.sh                      # all tests
tests/run.sh t_cli t_astrometry   # selected tests (the .sh suffix is optional)
KEEP=1 tests/run.sh t_gphot       # keep the working directories
```

Each test runs in its own empty directory, created by `mktemp -d` under `$TMPDIR` (default `/tmp`). The directory is removed after a pass and kept after a failure; the log gives its path.

### Output

Each check prints one line, so the log of a failed test shows everything that passed before the failure:

```text
Munipack 0.6.2+hg1653, (C) 1997-2026 F.Hroch ...
=== t_artificial
  ok: --seed reproducible
  ok: detector mean [ADU] = <value> (124.9 .. 125.1)
  ...
PASS t_artificial (<seconds> s)
...
=== 16 passed, 0 failed
```

A test stops at its first failed check, which is printed as `FAIL: <check>`. `run.sh` exits with status 0 only when every test passes.

## `make check` (upstream)

The automake test `src/test/check.sh` runs `find -f 5`, `aphot`, `astrometry` and `phcal` on the bundled frame `barnard-calibrated.fits`, with the catalogue `cone_barnard.fits` (Johnson V, `--area 0.3`).

- **Pass condition:** the calibration constant `CTPH` of the `PHOTOMETRY` extension is within twice its error `CTPHERR` of 11.5, with 0 < `CTPHERR` < 3.
- **Change from upstream:** upstream required 11 < `CTPH` < 12. That window is narrower than the error on this frame (about 1.3), so legitimate changes in the star selection made it fail.

The Docker build runs this test after `make install`. A failure stops the build and prints `test/check.sh.log`. In a configured and built source tree (see the `Dockerfile` for the `bootstrap`, `autoreconf` and `configure` steps), run `make check` in `src/`. The log is `src/test/check.sh.log`.

## The helper library `tests/lib.sh`

`run.sh` sources `lib.sh` before each test; `lib.sh` sets `set -e -o pipefail`.

| Helper | What it does |
|---|---|
| `fail MSG` | Prints `FAIL: MSG` to stderr and stops the test. |
| `ok MSG` | Prints `  ok: MSG`. |
| `check NAME VALUE LO HI` | Passes when VALUE is a number and LO ≤ VALUE ≤ HI, otherwise fails. A missing value (`''`) or `nan` fails. |
| `expect_fail NAME CMD...` | Passes when CMD exits with status 1: an error, not a crash (139) or a signal. |
| `key FILE[EXT] KEY` | The value of a header keyword (`munipack fits -K`), or `''` when the keyword is missing. |
| `nrows FILE[EXT]` | The number of rows of a table (`NAXIS2`). |
| `colnum FILE[EXT] NAME` | The 1-based index of a table column, or `''` when the column is missing. |
| `column FILE[EXT] NAME...` | The named columns of a table as whitespace-separated text. |
| `imstats FILE [x1:x2,y1:y2]` | Prints `mean sd n` of the image pixels, optionally of a section. |
| `imsum FILE` | The MD5 sum of the pixel values, used to compare images bit for bit. |
| `mkcat FILE.fits N SEED [RA DEC RADIUS MAGMIN MAGMAX MINSEP]` | Writes a synthetic catalogue (details below). |
| `catfits FILE.fits < lines` | Writes a catalogue in the format of `mkcat` from `RA DEC V B` lines on stdin. |
| `tabfits FILE.fits "NAME FORM UNIT ..." < lines` | Writes a binary table (extension `CATALOGUE`) with any columns. FORM is `D` or `E`; UNIT `-` means none. |
| `match CAT[EXT] "CATCOLS" TAB[EXT] "COLS" TOL` | Joins two tables by position. The first two columns of each are RA and Dec. Each table row is paired with its nearest catalogue star within TOL degrees, and the output is the remaining catalogue columns followed by the remaining table columns. The tests use TOL = 3e-4° (1.1″). |
| `stats < numbers` | Prints `median robust-sd n`, where the robust sd is 1.4826 × MAD. |
| `CATCOLS` | The column options for a `mkcat` catalogue: `--col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA --col-pm-dec pmDE --col-mag Vmag`. |

Details of `mkcat`:

- **Positions:** N stars at random positions within RADIUS degrees of RA, Dec, with no two stars closer than MINSEP degrees. The defaults are RA 180, Dec 0, RADIUS 0.08, MINSEP 0.004.
- **Magnitudes:** V is uniform from MAGMIN to MAGMAX (default 10 to 15), and B − V is uniform in 0 … 1.5.
- **Columns:** `RAJ2000`, `DEJ2000` [deg], `pmRA`, `pmDE` [mas/yr, zero], `Vmag`, `e_Vmag`, `Bmag`, `e_Bmag`. The errors are 0.001 mag, because `artificial` renders the magnitudes exactly.
- **Random numbers:** a Park–Miller generator inside awk, so the catalogue is the same with gawk and mawk (mawk ignores `srand(seed)`).
- **Failure:** the helper fails if it cannot place N stars.

`catfits` and `tabfits` write a FITS text dump and convert it with `munipack fits --restore`.

### Writing a test

A test is a bash fragment `tests/t_<area>.sh` that is sourced after `lib.sh`, so it needs no shebang. It starts with a comment that says what it simulates and which `BUGS.md` entries it covers. Commands that must succeed run under `set -e`, with their output sent to `/dev/null`. Every quantity is checked with `check`, `expect_fail` or an explicit `fail`. `run.sh` picks up new `t_*.sh` files automatically.

Conventions used below:

- `--scale` is in pixels per degree, so `--scale 3600` is 1″/px.
- `--hwhm` of `artificial --psf GAUSS` is the Gaussian sigma.
- Astrometric errors are given in pixels:
  - **CRVAL1** is the offset in RA times cos δ.
  - **CD** is the largest difference of a CD matrix element, expressed as the position error it causes 300 px from the reference pixel.

## The tests

| Test | Area | `BUGS.md` entries |
|---|---|---|
| `t_cli` | Front end: help, exit status, rejected options | A14, A16, A17, A24, P45 |
| `t_artificial` | Simulator: seeds, detector noise, random fields, `--mproc` | A2, A13, A19, A23, A25, A33, A34, §2.1 |
| `t_calib` | Master bias, master dark, `phcorr` | P31, P39, P41 |
| `t_pipeline` | `find` → `aphot` → `astrometry` → `phcal`, and `psfphot` → `phcal` | P8, P18, P29, P34, P35, P42, P46, S7 |
| `t_astrometry` | Solution from scratch on rotated frames at Dec +60 | S1–S3, S7, A1, A18 |
| `t_astrometry_wcs` | WCS-guided match, wrong sky, few stars, proper motions | S16, S21–S26 |
| `t_astrometry_roll` | Pointing match with unknown roll, threads, large cones | S29–S31 |
| `t_kombine` | Stacking: flux, positions, exposure time | S9–S14, S20 |
| `t_crowding` | Apertures beside brighter neighbours | P25, P29, P46 |
| `t_phcal_filters` | Two filters in any order | P14, P43 |
| `t_colterm` | Colour term of one filter, `phfotran --model mag`, `--tratab` | P47, P48 |
| `t_colsdss` | Colour term with real SDSS colours | P47, P55 |
| `t_colour2` | Two colour terms with real SDSS colours | P56 |
| `t_gphot` | Growth-curve photometry and its calibration | P6, P24, P38, P49 |
| `t_syserr` | `phcal --syserr auto` with few calibration stars | P54, P55 |
| `t_psfempir` | Empirical PSF with radial variation | P50–P53 |

### `t_cli`: front end

**Simulates:** nothing; the test runs the front end only.

**Checks:**
- `munipack --version`, `munipack --help`, and `munipack <action> --help` succeed for every action: `bias dark flat phcorr kombine astrometry find aphot gphot psfphot phcal artificial phfotran timeseries fits cross`.
- Exit status exactly 1 for each of these:
  - an unknown action;
  - a missing input file (`find`);
  - `artificial --psf FOO`;
  - `artificial --spread FOO`;
  - `artificial --mproc 0`;
  - `artificial --rcen` without `--dcen`.
- A rejected `artificial` command writes no output file.
- `munipack phfotran` without input fails with an error about its input, not "unknown action". This shows that `phfotran` is dispatched.

**Covers:**
- A24: errors exited with 139.
- A14: unknown PSF or spread names were silently replaced.
- A16: paired options were silently dropped.
- A17: `--mproc ≤ 0` was accepted.
- P45: `phfotran` was not dispatched.

### `t_artificial`: frame simulator

**Simulates:** 200 × 200 px frames of an empty field (one star of 30 mag), plus random star fields without a catalogue.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| `--seed 7` twice, `--seed 8` | same seed gives a bit-identical frame; another seed gives a different frame | exact |
| Dark 0.5 e⁻/s/px, 100 s, read noise 5 e⁻, gain 2, bias 100 ADU | mean = 100 + 50/2 = 125 ADU | 124.9 … 125.1 |
| | sd = √(50 + 25)/2 = 4.330 ADU | 4.24 … 4.42 |
| | `GAIN` = 2, `RDNOISE` = 5 | exact |
| The same with `--disable-noise` | mean 125 ADU | 124.999 … 125.001 |
| | sd | 0 … 0.01 |
| `--sky-mag 99` alone | mean (no sky, no dark) | 0 … 0.001 |
| `-B 16` | mean (rounded; truncation would give 124.5) | 124.9 … 125.1 |
| | `BITPIX` = 16 | exact |
| Random field, default 666² px | stars in `CATALOGUE` | 1 … 20 |
| Random field, 1332² px (4× the area) | stars in `CATALOGUE` | 15 … 60 |
| Two frames, `--mproc 2` against one process | frames bit-identical | exact |

**Covers:**
- A2: reproducibility by `--seed`, independent of `--mproc`.
- A13, A19 and §2.1: the detector-noise model.
- A33: `--sky-mag 99` left a sky.
- A34: integer frames were truncated.
- A23: the star density did not depend on the frame area.
- A25: `--mproc`.

### `t_calib`: calibration frames

**Simulates:** 150 × 150 px frames of a low-count, space-UV-like detector: gain 1, read noise 0.4 e⁻, pedestal 100 ADU. The test makes 10 bias frames and 10 dark frames of 100 s at 0.0013 e⁻/s/px (0.13 e⁻/px), plus a science frame with a 21 mag sky and dark.

**Checks:**

| Check | Tolerance |
|---|---|
| Master bias (`bias`), mean | 99.995 … 100.005 ADU |
| Master dark (`dark -bias`): `DARKMEAN` mentions `CLIP` | exact |
| Master dark, mean (true 0.13 e⁻/px) | 0.127 … 0.133 |
| Master dark, `EXPTIME` | 100 |
| `phcorr -bias -dark`: level minus the noiseless sky (3638.5 e⁻/px) | −1 … +1 e⁻/px |
| `phcorr`: `GAIN` = 1, no `BIASLEV` (the output is in electrons) | exact |

**Covers:**
- P31: the robust dark mean was about 16% low at low counts.
- P39: the noise model of `phcorr`.
- P41: the bias of a single dark frame was not subtracted.

### `t_pipeline`: the main route

**Simulates:** 150 stars with V = 10–16 at RA 180, Dec +30, on a 640² px frame at 1″/px. The PSF is Gaussian (σ = 1.5 px); the detector has gain 2, read noise 5 e⁻, bias 100 ADU and dark 0.01 e⁻/s/px.

**Checks:**

| Step | Check | Tolerance |
|---|---|---|
| `find -f 3.5` | stars found | 145 … 155 |
| `aphot` | `HWHM` of `APERPHOT` (true 1.1774 σ = 1.766 px, pixelised) | 1.65 … 1.95 px |
| `astrometry --remove`, then `astrometry -c` | `CRVAL1` removed by `--remove` | none left |
| | `CRPIX1` | 320 |
| | `CRVAL1`, `CRVAL2` offsets | ±0.05 px |
| | scale | 3599.5 … 3600.5 pix/deg |
| | rotation `CD1_2/CD2_2` | ±3e-4 |
| | `CD1_1` < 0 (east left) | exact |
| `phcal` (Johnson V) of `aphot` | stars matched; valid `MAG` | 145 … 150 each |
| | median `MAG − V` | ±0.003 mag |
| | robust sd of `MAG − V` | ≤ 0.005 mag |
| | robust sd of `(MAG − V)/MAGERR` | 0.7 … 1.4 |
| `psfphot` on the frame with the original WCS | stars in `PSFPHOT` | 145 … 155 |
| `phcal --photometry psf` | `PHOTOTYP` mentions `PSF`; the same calibration checks as for `aphot` | as above |

**Covers:**
- P8: `HWHM` stored σ.
- S7: RA residuals without cos δ.
- P29 and P46: isolated stars kept their apertures.
- P34: aperture errors ignored the gain.
- P18, P35, P42: `psfphot` speed, variance and faint-star errors.

### `t_astrometry`: rotated frames at Dec +60

**Simulates:** 80 stars within 0.08° of RA 30, Dec +60, where cos δ = 0.5. Frames are 600² px at 1″/px with a Gaussian PSF (σ = 1.2 px), rotated by `--angle 30` and `--angle -120`. The WCS of `artificial` is removed, and `astrometry -c` solves each frame from scratch.

**Checks:**
- `CRVAL1` and `CRVAL2` offsets from the truth: ±0.05 px.
- CD matrix error against the true CD: ≤ 0.05 px at 300 px.
- `astrometry --mproc 2` gives keywords identical to the one-process solution: `CRVAL1/2` and `CD1_1` to `CD2_2`.

**Covers:**
- S1–S3: rotation, reflection and the CD/PC matrix.
- S7: RA residuals with cos δ.
- A1, A18: `--scale` is in pix/deg.

### `t_astrometry_wcs`: matching by the header WCS

**Simulates:** 150 stars within 0.12° of RA 30, Dec +60, on 600² px frames at 1″/px rotated by 20°. Unless stated otherwise, every solution is checked against the true WCS with the tolerances CRVAL ±0.05 px and CD ≤ 0.05 px at 300 px.

**Checks:**
1. **Perturbed header.** The header's `CRPIX` is moved by (+45, −60) px (75 px), the CD matrix is rotated by 2.5° and scaled by 2.5%. The frame is solved with `ASTRMETH = 'WCS'`.
2. **Blind match.** The same frame with `--blind` is solved with `ASTRMETH = 'BLIND'`.
3. **Wrong sky.** The catalogue holds different stars at the same place. The match by the header and `--blind --scale-free` both exit with status 1, and no WCS is saved (`ASTRMETH` missing).
4. **Small frame.** A 200² px frame (about 15 stars) with a 30-star catalogue of its footprint and the WCS removed is solved, with 5–40 stars used. It is also solved with `--scale 3600`.
5. **Few stars.**
   - `--minmatch 3` is rejected with exit status 1.
   - A run of three frames whose middle frame has only four stars exits non-zero with "too few stars".
   - The other two frames are still calibrated (`ASTRMETH = 'WCS'`).
6. **Proper motions.**
   - Proper motions of ±400 mas/yr, with a catalogue of epoch 2000 (`JEPOCH` in the table header) and a frame dated 2030-01-01. The frame is solved within the tolerances.
   - With a catalogue that has proper motions but no epoch, the frame is still solved, with the warning "proper motions of … ignored" on stderr and "proper motions ignored" in the header.
7. **Catalogue order.** The catalogue in reversed order gives a solution within the tolerances.

**Covers:**
- S21: the header WCS was ignored.
- S22: the scale prior.
- S16, S23: the significance test against false solutions.
- S24: few stars and multi-frame runs.
- S25: proper motions with cos δ and a warning when the epoch is missing.
- S26: catalogue order.

### `t_astrometry_roll`: pointing known, roll unknown

**Simulates:** the same field and frames as `t_astrometry_wcs`, with the same solution tolerances. The header WCS is perturbed by a shift, a rotation, a scale factor and an optional mirror; RA is shifted via the old CD matrix.

**Checks:**
1. **Rotated header.** The header is rotated by 137° (beyond `--wcs-rot`), moved by (30, −40) px and scaled by 1.015. The frame is solved with `ASTRMETH = 'POINTING'`.
2. **`--roll-known`.** On the same header, `--roll-known` trusts the wrong rotation, falls back to the blind match, and is solved with `ASTRMETH = 'BLIND'`.
3. **Pointing without a header WCS.** The WCS is removed and the frame is solved with `--rcen`, `--dcen` (about 30″ off) and `--scale 3564` (1% off): `POINTING`.
4. **Mirrored header.** The header is mirrored, rotated by 75° and scaled by 0.99. `--reflex-unknown --roll-unknown` solves it: `POINTING`.
5. **Units of `--wcs-offset`.**
   - The header is moved by 250 px and rotated by −100°, and `--wcs-offset 5arcmin` solves it: `POINTING`.
   - `--wcs-offset 3parsec` exits with status 1.
   - `--roll-known --roll-unknown` together exit with status 1.
6. **Wrong sky.** With the widest search (`--roll-unknown --reflex-unknown --wcs-offset 50%`), the wrong sky exits with status 1 and no WCS is saved.
7. **Large cone.** A cone of 0.3° with 900 stars (3.4× the frame), with the pointing given by `--rcen`, `--dcen` and `--scale 3600`, is solved within the tolerances, using 60–250 stars.
8. **Threads.** The blind match (`--blind --scale 3600`) of the cone frame without a pointing, run with `OMP_NUM_THREADS=1` and with `OMP_NUM_THREADS=4`, gives the same `ASTRMETH` and identical WCS keywords. The check does not require the match to succeed.

**Covers:**
- S29: matching by pointing and scale with unknown roll.
- S30: the OpenMP blind match is independent of the number of threads.
- S31: the catalogue is selected near the frame.

### `t_kombine`: stacking

**Simulates:** 40 isolated stars (at least 0.008° = 29″ apart) within 0.06° of RA 100, Dec +10, on 500² px frames at 1″/px (Gaussian σ = 1.3 px, gain 1, read noise 3 e⁻, 30 s). There is one reference frame and four frames shifted in RA by about 5–16 px. `kombine -i bilinear` stacks the four frames onto the reference grid.

**Checks:**

| Check | Tolerance |
|---|---|
| `EXPTIME` of the stack | 120 s |
| Stars found in both the reference and the stack (within 1 px) | 36 … 40 |
| Median `APCOUNT8` ratio, stack/reference | 3.96 … 4.04 |
| Median x and y offset of the stars | ±0.03 px |

**Covers:**
- S20: the stack is the sum.
- S14: the flux is preserved.
- S9–S13: the resampling grid and the positions given by the WCS.

### `t_crowding`: apertures beside brighter neighbours

**Simulates:** designed pairs of stars on a 400² px frame: pairs 12–20 px apart, and a 14 mag star 7 px from a 9 mag one. Isolated stars give the frame's mean profile. The frame is made twice: with a Gaussian PSF (σ = 1.5 px) and with a Moffat PSF (`--hwhm 2.7 --beta 1.5`, heavy wings).

**Checks (both PSFs):**
- The 14 mag star 7 px from the 9 mag one has an invalid `APCOUNT7` (5.1 px aperture), or is not detected at all (`find` rejects detections in the wings, P25).
- Three isolated stars and the 9 mag member of that pair have a valid `APCOUNT7`.

**Checks (Gaussian only):** the fainter member of each 12–20 px pair (V = 11.5, 12.5, 11.0) has `APCOUNT7` equal to that of an isolated 11 mag star, scaled by magnitude, within 0.995 … 1.005.

**Covers:**
- P29: a neighbour's light stayed in the apertures.
- P46: a star's own light in a neighbour's growth curve invalidated its apertures.
- P25: secondary detections.

### `t_phcal_filters`: two filters

**Simulates:** 80 stars (B − V = 0 … 1.5) within 0.07° of RA 200, Dec −20. B and V frames are rendered from the catalogue B and V, 540² px, Gaussian σ = 1.4 px.

**Checks:** `phcal` with `-f B,V` and with `-f V,B` (columns and files in the same order as the filters). For each filter:

| Check | Tolerance |
|---|---|
| `FILTER` of the output | exact |
| `PHOTOMETRY` tables of the two orders | identical |
| Stars matched with a valid `MAG` | 76 … 80 |
| Median `MAG − catalogue` | ±0.003 mag |
| Robust sd | ≤ 0.004 mag |
| Slope of the residuals against B − V | ±0.003 mag/mag |

**Covers:**
- P14: the order of the filters.
- P43: a multi-filter calibration introduced a colour term.

### `t_colterm`: the colour term of one filter

**Simulates:** an instrument whose bands differ from the catalogue's: B_i = B − 0.10 (B − V) and V_i = V + 0.04 (B − V). It uses 80 stars and B and V frames as in `t_phcal_filters`.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| `phcal -f B --color-term Bmag,Vmag,e_Bmag,e_Vmag` | `COLTK1` (true 0.10) | 0.095 … 0.105 |
| | `COLORSRC` column present | exact |
| | stars | 76 … 80 |
| | slope of the residuals against B − V | ±0.003 |
| `phcal` without `--color-term` | no `COLTK1` keyword | exact |
| `--color-term` with `-r` | exit status 1 | exact |
| `phfotran --model mag -f B,V` | `TRAMODEL` = `MAG` | exact |
| | K1 of B (0.10) | 0.095 … 0.105 |
| | K1 of V (−0.04) | −0.045 … −0.035 |
| `phcal --tratab` with this table, B and V | robust sd | ≤ 0.004 mag |
| | colour slope | ±0.003 |

**Covers:**
- P47: the new `--color-term` option.
- P48: the magnitude model of `phfotran`.

### `t_colsdss`: colour term with real colours

**Simulates:** 93 stars of SDSS DR16 (VizieR V/154/sdss16; primary, clean, class star), embedded in the test. They lie within 0.27 × 0.27° around RA 180, Dec 0, with 14 < g < 20 and errors below 0.05 mag. A 1024² px frame in Gunn g is rendered from the catalogue g. It is calibrated against the catalogue r with the colour g − r, so the true term is m_cat − m_inst = ZP + k1 (g − r) with k1 = −1 exactly.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| `--color-term gpmag,rpmag,e_gpmag,e_rpmag` | `COLTK1` | −1.02 … −0.98 |
| | `COLTK1E` | 0.001 … 0.01 |
| | `(COLTK1 + 1)/COLTK1E` | ±3 |
| | stars matched | 85 … 93 |
| | slope of `MAG − r` against g − r | ±0.015 |
| | robust sd of `MAG − r` | ≤ 0.015 mag |
| | `PHSYSERR` | ≤ 0.01 mag |
| `--color-order 2` | `COLTK1` | −1.02 … −0.98 |
| | `COLTK2` | ±0.02 |
| No colour term | slope against g − r | 0.95 … 1.05 |
| | `PHSYSERR` (the colours as scatter) | 0.3 … 0.5 mag |

**Covers:**
- P47: the colour term.
- P55: the errors of the colour coefficients when the systematic error is estimated.

### `t_colour2`: two colour terms

**Simulates:** 61 SDSS DR16 stars embedded in the test, selected as in `t_colsdss` but with 14.5 < g < 19.5 and e_u < 0.08. The colours u − g and g − r are correlated by 0.96. A 1024² px frame in Gunn g is rendered from X = g + 0.3 (u − g) − 0.2 (g − r) and calibrated against g with both colours, so the true coefficients are k1 = −0.3 and k2 = +0.2. Three stars have no u in the calibration catalogue; their u − g is predicted from g − r.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| `--color-term upmag,gpmag --color-term gpmag,rpmag` | `COLTNTRM` | 2 |
| | `COLTERM2` | `gpmag-rpmag` |
| | `COLTK1` (−0.3) | −0.33 … −0.27 |
| | `COLT2K1` (+0.2) | 0.15 … 0.25 |
| | pulls `(COLTK1 + 0.3)/COLTK1E`, `(COLT2K1 − 0.2)/COLT2K1E` | ±3 |
| | correlation `COLTC12` | −1 … −0.5 |
| | `COLOR2ERR` column present | exact |
| | stars matched | 55 … 61 |
| | robust sd of `MAG − g` (stars with both colours, `COLORSRC` 3) | ≤ 0.012 mag |
| | slope against u − g; against g − r | ±0.01; ±0.02 |
| | stars without u (`COLORSRC` 2) | 2 … 3 |
| | their largest \|`MAG − g − median`\|/`MAGERR` | ≤ 3 |
| The same with the colour errors (4-element `--color-term`) | both pulls | ±3 |
| One `--color-term upmag,gpmag` | no `COLTNTRM` | exact |
| | robust sd | 2× the two-colour sd … 0.1 mag |
| Dependent colours (`upmag,gpmag` and `gpmag,upmag`) | exit status 1 | exact |
| Four colour terms | exit status 1 | exact |
| `--color-order 1,2,1` for two terms | exit status 1 | exact |

**Covers:** P56, several colour terms.

### `t_gphot`: growth-curve photometry

**Simulates:** 200 stars with V = 10–16.5 within 0.1° of RA 180, Dec +30, on an 800² px frame (Gaussian σ = 1.5 px, gain 2, read noise 5 e⁻, dark 0.01 e⁻/s/px, sky 21 mag). `aphot` runs with 11 apertures from 1 to 18 px and a sky ring of 20–30 px. The truth is the `CATALOGUE` extension of `artificial`, where `COUNT/GAIN` is the total flux in ADU.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| `gphot` | `GROWMETH` = `EE` | exact |
| | `EENSTAR` (stars of the growth curve) | 30 … 200 |
| | `HWHM` of `GROWPHOT` (true 1.766 px) | 1.6 … 2.1 px |
| | valid `GCOUNT` and `GCOUNTERR` | 97 … 100% |
| | stars joined with the truth (1 px) | 180 … 200 |
| | median `GCOUNT/truth − 1` at S/N > 50 | ±0.003 |
| | median pull `(GCOUNT − truth)/GCOUNTERR` | ±0.3 |
| | robust sd of the pull | 0.8 … 1.25 |
| | stars with S/N > 300 | 15 … 200 |
| | robust sd of their pull | 0.7 … 1.3 |
| | median `(SKY − SKY of APERPHOT)/SKYERR of APERPHOT` | ±0.5 |
| | stars with `SKYERR` larger than that of `APERPHOT` | 0 |
| `phcal --photometry grow` | `PHOTOTYP` mentions `GPHOT` | exact |
| | stars calibrated | 180 … 200 |
| | median `MAG − V` (`MAGERR` < 0.01) | ±0.003 mag |
| | robust sd `MAG − V` (`MAGERR` < 0.01) | ≤ 0.006 mag |
| | robust sd of `(MAG − V)/MAGERR` | 0.7 … 1.3 |
| `gphot --legacy` | `GROWMETH` = `LEGACY` | exact |
| `gphot -th 0.2` without `--legacy` | exit status 1 | exact |
| `gphot` on a frame of two stars | exit status 1 (not a crash) | exact |

**Covers:**
- P24 and P49: GCOUNT is the enclosed energy of the frame's growth curve.
- The errors of the growth curve are included in `GCOUNTERR`.
- P6 and P38: the corrected sky and its error.

### `t_syserr`: systematic error from few stars

**Simulates:** 30 stars (V = 11–14) on a 300² px frame rendered from the catalogue. The calibration catalogues scatter around the truth:
- `few.fits`: 6 stars, by about 0.05 mag.
- `many.fits`: all 30 stars, by about 0.03 mag.

**Checks:**

| Case | Check | Tolerance |
|---|---|---|
| 6 stars, `--syserr auto` (default) | `PHSYSNST` | 6 |
| | `PHSYSERR` | 0.02 … 0.10 mag |
| | `PHSYSERE` (upper error) | 0.01 … 0.10 mag |
| | warning "estimated from 6 stars only" on stderr | once |
| | `COMMENT`s "68% interval" and "Warning: the systematic" | once each |
| | `1.0857 CTPHERR/CTPH × √6 / PHSYSERR` (at least t(5) = 1.11; it was 1.0) | 1.06 … 3 |
| `--syserr 0.05` | `PHSYSERR` | 0.05 |
| | no `PHSYSERE`, no `PHSYSNST`, no warning | exact |
| 30 stars | `PHSYSNST` | 28 … 30 |
| | `PHSYSERR` | 0.02 … 0.04 mag |
| | `PHSYSERE` | 0.002 … 0.01 mag |
| | no warning | exact |

**Covers:**
- P54: the systematic error estimated from few stars, with its upper error and a warning.
- P55: `CTPHERR` includes the uncertainty of the estimated systematic error.

### `t_psfempir`: empirical PSF

**Simulates:** 150 stars with V = 10–16 at RA 180, Dec +30, on a 640² px frame with a wide-winged Moffat PSF (`--hwhm 2.5 --beta 2.5`) and a ground-like detector. The test then runs `psfphot --psf EMPIRICAL --psf-var-radial 1`, followed by `phcal --photometry psf`.

**Checks:**

| Check | Tolerance |
|---|---|
| Stars in `PSFPHOT` | 145 … 155 |
| `PSFVARR` | 1 |
| `PSFRAD/FWHM` (99.5% of the flux of β = 2.5 is within 5.1 FWHM; it was 4) | 4.5 … 5.8 |
| `PSFNORMN` (stars of the normalisation) | 10 … 100 |
| `PSFSMOOT` (the regularisation) | 0.001 … 10 |
| `PSFMODERR` column present and ≤ `PSFCOUNTERR` for every star | exact |
| `phcal`: stars matched | 140 … 150 |
| Median `MAG − V` | ±0.003 mag |
| Robust sd `MAG − V` | ≤ 0.006 mag |
| Robust sd of `(MAG − V)/MAGERR` | 0.7 … 1.4 |

**Covers:**
- P50: the table radius by the growth curve and normalisation to an aperture.
- P51: radial variation.
- P52: `PSFMODERR` in the errors.
- P53: cross-validated regularisation.
