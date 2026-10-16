# astrometry — astrometric calibration

`astrometry` computes the astrometric calibration (WCS) of a frame and writes
it to the primary header. It works in two steps. First it **matches** the
detected stars to a reference catalogue or to an already calibrated frame,
which pairs each frame star with its catalogue star. Then it **fits** the
projection (gnomonic `TAN`, or `TPV` with polynomial distortion) to those
pairs by a robust fit. When the frame header already holds an approximate
WCS, or when you give the pointing (`--rcen`, `--dcen`) and the scale, the
match is guided by it. Otherwise the match is blind: it searches sequences
of star triangles. Every solution is tested against chance coincidences
before it is saved. Frames are independent: when one frame fails, the rest
of the run is still processed, and the exit status is 1.

## Synopsis

```bash
munipack astrometry [options] file(s)
munipack astrometry [options] 'input.fits output.fits'
```

By default each file is updated in place. To write the result to a new file,
give the input and the output name as one argument separated by a space
(quoted, as above), or use `-t`/`-O`. A list of files (one name, or an
input and an output name, per line) can be read from a file with `@list.txt`
or from standard input with `-` as the file name.

In a container:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
    astrometry -c cone.fits frame.fits
```

**Prerequisites.** The matching modes (`match`, `cross`, `sequence`) need:

- the `FIND` extension from [`find`](find.md) (required);
- the `APERPHOT` extension from [`aphot`](aphot.md) (strongly recommended).
  Its first aperture ranks the stars by brightness, feeds the flux check of
  the match and gives the weights of the final fit. Without it all stars get
  the same flux: they stay in the order of `FIND` and the fit is unweighted.

They also need a reference: a catalogue (`-c`, e.g. from [`cone`](catalogues.md))
or a calibrated frame (`-r`, `-R`). This build has no default catalogue.
Without a reference the run fails with "too few stars … in the reference".
The `manual` mode needs none of these.

## Modes

| Mode | What it does | ASTRMETH |
|---|---|---|
| `match` (default) | Guided match by the WCS in the header or by the pointing, falling back to the blind match; then the fit. | `WCS`, `POINTING` or `BLIND` |
| `cross` | Cross-matches the stars by an approximate calibration that is known within the tolerances of the guided match. The calibration comes from the header WCS or from `--scale`, `--angle`, `--reflex` (and `--rcen`, `--dcen`). The offset is found by the same vote as in the guided match. If the vote fails, the nearest stars within 5 × `--sig` are used. Add `--roll-unknown` to allow any rotation. | `CROSS` |
| `sequence` | Uses a correspondence that you give: `--seq1` lists catalogue (or reference) indices and `--seq2` the matching frame indices. Both lists must have the same length. | `SEQUENCE` |
| `manual` | Writes the WCS from `-p`, `--rcen`, `--dcen`, `--xcen`, `--ycen`, `--scale`, `--angle`, `--reflex` with no catalogue, no star table and no checks. | — |

Mode names are case-insensitive.

## Options

### Mode and reference

| Option | Argument | Default | Description |
|---|---|---|---|
| `-m`, `--mode` | `match`, `cross`, `sequence`, `manual` | `match` | Calibration mode (see above). |
| `-c`, `--cat` | FITS table | — | Reference catalogue. FITSIO filters can be used, e.g. `'cone.fits[1][$f.mag$>12]'`. |
| `-r`, `--ref` | FITS file | — | Already calibrated frame. Its stars are the reference. |
| `-R`, `--rel` | FITS file | — | Reference frame for relative astrometry: no projection, coordinates in pixels. |
| `--vocat` | XML file | built-in `VOcat_conf.xml` | VO catalogue configuration. It maps the catalogue to its default column labels. |
| `--col-ra` | label | by the catalogue (`RAJ2000`) | Right Ascension column [deg]. |
| `--col-dec` | label | by the catalogue (`DEJ2000`) | Declination column [deg]. |
| `--col-pm-ra` | label | by the catalogue (`pmRA`) | Proper-motion column in RA, μ<sub>α</sub> cos δ. RA moves by pmRA·Δt / cos δ. |
| `--col-pm-dec` | label | by the catalogue (`pmDE`) | Proper-motion column in Dec. |
| `--col-epoch-ra` | label | by the catalogue | Reference epoch of RA, per star. |
| `--col-epoch-dec` | label | by the catalogue | Reference epoch of Dec, per star. |
| `--col-mag` | label | by the catalogue (`f.mag` for UCAC) | Magnitude-like column. The catalogue is sorted by it. |

Give only one of `-c`, `-r` and `-R`. Proper motions are applied from the
epoch of the catalogue (the table key `JEPOCH` or the epoch columns) to the
`DATE-OBS` of the frame. If the catalogue has proper motions but no epoch,
they are ignored with a warning on stderr and a `COMMENT` in the header.

### Projection, pointing and manual parameters

| Option | Argument | Default | Description |
|---|---|---|---|
| `-p`, `--projection` | `TAN`, `TPV` | `TAN` | `TAN`: gnomonic. `TPV`: gnomonic with third-order polynomial distortion (`PV1_0…PV2_11`). |
| `--xcen` | real [pix] | width/2 | Reference pixel x (CRPIX1). Give it together with `--ycen`; one of them alone is ignored. |
| `--ycen` | real [pix] | height/2 | Reference pixel y (CRPIX2). |
| `--rcen` | real [deg] | — | RA of the reference pixel (CRVAL1). In `match` mode, together with a scale, it is the **pointing** for the guided match. Give it together with `--dcen`; one of them alone is ignored. |
| `--dcen` | real [deg] | — | Dec of the reference pixel (CRVAL2); the pointing in `match` mode. |
| `--scale` | real [pix/deg] | — | Plate scale **in pixels per degree** (1″/px = 3600). In `match` mode it is the scale prior. |
| `--angle` | real [deg] | — | Rotation, clockwise positive. In `match` mode it is the known roll for the guided match. |
| `--reflex` | — | off | The frame is mirrored. |
| `--radesys` | string | `ICRS` | Reference system written to `RADESYS`. |
| `--epoch` | real [yr] | from the catalogue | Reference epoch written to `JEPOCH`. |

### Guided match and scale prior

| Option | Argument | Default | Description |
|---|---|---|---|
| `--blind` | — | off | Ignore the header WCS and the pointing: no guided match and no scale prior from the header. `--scale`, if given, is still used as the prior. |
| `--scale-tol` | real | 0.05 | Relative tolerance of the scale prior, which comes from `--scale` or the header WCS. 0 or less: any scale. |
| `--scale-free` | — | off | Match with any scale. Same as `--scale-tol 0`. Use it when the header scale is wrong by more than the tolerance. |
| `--wcs-offset` | number[unit] | 10% of the frame, ≥ 100 pix | Tolerance of the offset of the header WCS or of `--rcen`, `--dcen`. A plain number is pixels; units `pix`, `arcsec` (`"`), `arcmin` (`'`), `deg`, `%` (of the larger side) are accepted, e.g. `90arcsec`, `3arcmin`, `0.1deg`, `20%`. The vote uses at most the frame diagonal. |
| `--wcs-rot` | real [deg] | 3 | Tolerance of the rotation of the WCS or of `--angle`. 180 or more means any rotation, as `--roll-unknown`. |
| `--wcs-scale` | real | 0.03 | Relative tolerance of the scale used by the guided match. |
| `--roll-unknown` | — | off | Search any rotation at once and skip the pass within `--wcs-rot`. |
| `--roll-known` | — | off | The rotation is known within `--wcs-rot`: no pass over any rotation. If that fails, the blind match follows. Cannot be combined with `--roll-unknown`. |
| `--reflex-unknown` | — | off | The guided match tries both reflections, which doubles the trials. |

### Blind match and acceptance tests

| Option | Argument | Default | Description |
|---|---|---|---|
| `--sig` | real [pix] | 1 | Error of the frame centroids. It sets the match tolerances, the vote bins (4 × `--sig`) and the cross-match limit (5 × `--sig`). |
| `--sigcat` | real [deg] | 1″ (2.8e-4) | Error of the catalogue or reference coordinates. |
| `--flxmin` | real | 0.63 (−0.5 mag) | Lower limit of the flux ratio of matched stars. |
| `--flxmax` | real | 1.58 (+0.5 mag) | Upper limit of the flux ratio of matched stars. Widen the range for a filter that differs from the catalogue band. |
| `--minmatch` | integer ≥ 5 | 5 | Minimum length of the match sequence. Values below 5 are an error. Raise it to 7 or more on crowded fields. |
| `--maxmatch` | integer | 42 | Stars per chunk. Chunks of the frame and of the catalogue are matched against each other. 30–50 suits crowded fields. |
| `--maxstars` | integer | 250 | Number of the brightest stars used from the frame and from the catalogue. When the pointing is known, the catalogue stars are taken near the frame. |
| `--luckymatch` | integer | 3 | The match is accepted as reliable at `--minmatch` + this many stars. |
| `--disable-lucky-match` | — | off | Stop at the first successful match. Same as `--luckymatch 0`. |
| `--disable-flux-check` | — | off | Do not use fluxes during matching. Mismatches become more likely on crowded fields. |
| `--disable-rms-check` | — | off | Skip the RMS check and the significance test. Use with care. |

### Sequence mode

| Option | Argument | Default | Description |
|---|---|---|---|
| `--seq1` | `i,j,k,…` | — | Indices of the reference (catalogue) stars. |
| `--seq2` | `i,j,k,…` | — | Indices of the corresponding frame stars. Same length as `--seq1`. |

### Fit and saving

| Option | Argument | Default | Description |
|---|---|---|---|
| `--fit` | `robust`, `squares` | `robust` | Fitting method. |
| `--units` | `deg`, `arcmin`, `arcsec`, `mas`, `uas` | auto | Units of the RMS in the header comments. By default they follow the size of the RMS (pixels for `-R`). `pix`, listed by `--help`, is ignored (auto). |
| `--disable-save` | — | off | Compute the calibration without saving it. |
| `--remove` | — | off | Remove the calibration: the WCS keys, `PVi_j`, `RADESYS`, `JEPOCH`, `ASTRMETH` and the comment block. |
| `--report` | — | off | Add the `ASTROREPORT` extension (see Outputs). |
| `--show-defaults` | — | — | Print the default parameters and exit. |

### Output files

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Write the output files into this directory. |
| `-O` | — | off | Turn on advanced output filenames (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against input filenames (with `-O`). |
| `--mask` | string | — | Output filename built from back-references `\0`, `\1`, … (with `-O`). |
| `--format` | printf format | — | Format for numbered output files, e.g. `out_%02d.fits` (with `-O`). |

### Parallel processing and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-mproc`, `--mproc` | integer | 1 | Number of parallel processes; the frames are shared among them. The solutions do not depend on N. |
| `--Xproc` | string | — | Extra options passed to `mpirun`. |
| `--verbose` | — | off | Print the match, the vote, the significance, the number of threads and the fit. |
| `-h`, `--help` | — | — | Show the help and exit. |

## How the match works

**Order of the passes (`match` mode).** The guided match is used when the
frame has a WCS in its header (CD, or CDELT with PC), or when you give
`--rcen`, `--dcen` with a scale from `--scale` or the WCS. Unless `--blind`
is set, the passes are tried in this order:

1. **Known rotation** (when the WCS or `--angle` gives one; skipped with
   `--roll-unknown`). The catalogue is predicted on the frame and rotated and
   scaled on a grid within `--wcs-rot` and `--wcs-scale`. At each node of
   the grid, the position differences of the 60 brightest frame stars and the
   nearby catalogue stars vote for the offset, within `--wcs-offset`.
   `ASTRMETH = 'WCS'`.
2. **Any rotation** (skipped with `--roll-known`). Pairs of the 30 brightest
   frame stars and the catalogue stars of matching length do not depend on
   the rotation or the offset. Each such pair votes for a rotation and an
   offset. The vote above is then repeated around the eight best candidates.
   `ASTRMETH = 'POINTING'`. This pass is the right one when the pointing and
   the scale are known but the roll is not.
3. **Blind match**, with the scale prior. It is tried first on the
   catalogue stars near the frame (by the pointing), then on the brightest
   stars of the whole catalogue. Triangle sequences of at least `--minmatch`
   stars are built from chunks of `--maxmatch` stars. `ASTRMETH = 'BLIND'`.

**Scale prior.** When the scale is known from `--scale` or the header WCS,
the blind match skips star pairs whose length disagrees with it by more than
`--scale-tol` (5%). Solutions outside the tolerance are rejected. This makes
matching tens of times faster and removes false solutions with a wrong scale.
Without any scale the match is scale-free.

**Significance test.** Every solution is tested against chance. The test
counts the stars within 3 × the scatter of the fit, apart from the two that
fix the transformation. That count is compared with the Poisson number of
coincidences expected from the catalogue density in the overlap of the
catalogue and the frame. The probability is multiplied by the number of
trials of the search, and the solution is accepted below 10⁻³. Larger
searches (big `--wcs-offset`, `--reflex-unknown`, any rotation) count more
trials, so a frame of only 5–6 stars needs a tight `--wcs-offset`.
`--verbose` prints the numbers on the `# significance:` line.

**Five stars minimum.** The solution has four parameters (scale, rotation,
centre) and their errors, so at least five matched stars are required. A
frame or reference with fewer stars fails with "too few stars". Seven or
more stars is safer, and over 20 is ideal.

**Final fit.** The TAN fit of ten stars or more is weighted by the centroid
errors, estimated from the S/N in the first aperture and the PSF width, plus
an excess scatter. With `-p TPV` the fit goes in stages:

1. the gnomonic solution;
2. a second-order deformation, with all stars cross-matched again by it until
   the matched set is stable;
3. the third order, iterated the same way.

Distortions of about 15 px are recovered to about 0.1 px. With 25 matched
stars or fewer, TPV falls back to TAN with a warning. The TPV terms are fitted
without weights. Atmospheric refraction is not modelled.

**Threads.** The blind match runs the chunk pairs in OpenMP threads, and the
solution and significance are the same for any number of threads. The thread
count is `OMP_NUM_THREADS` (in a container pass it with
`docker run -e OMP_NUM_THREADS=4 …`). When it is unset, the count is the CPUs available
to the process (its affinity mask and CPU quota, e.g. `docker run --cpus 4`)
divided by `--mproc`. `--verbose` prints the count. Threads help most when no
solution exists, for example a wrong field without a scale prior.

## Choosing the catalogue cone

- **Pointing known** (header WCS, or `--rcen`, `--dcen`, `--scale`): the
  whole catalogue is read, and the brightest `--maxstars` stars within the
  frame half-diagonal + 10% + `--wcs-offset` are used. A cone larger than
  the frame is fine.
- **No pointing:** the brightest `--maxstars` stars of the whole cone are
  used. Make the cone about the field of view plus a 10% margin. For a square
  frame of side *L*, a radius of about 0.8 *L* works (the half-diagonal
  0.71 *L* plus 10%).
- **Magnitudes:** match the magnitude range to the stars that are well
  exposed and unsaturated on the frame. Use `cone --magmin/--magmax` or a
  FITSIO filter, e.g. `-c 'cone.fits[1][$f.mag$>12]'`. Sorting the cone is not
  needed: `astrometry` sorts it by `--col-mag`.
- **Crowded fields:** use `--minmatch 7 --maxmatch 44`, and on large
  telescopes also limit the bright end of the catalogue.

## Outputs

The calibration is written to the primary header and replaces any earlier
one. The image data are not changed.

| Key | Content |
|---|---|
| `CTYPE1`, `CTYPE2` | `RA---TAN`, `DEC--TAN` (or `…TPV`) |
| `CRVAL1`, `CRVAL2` | RA, Dec of the reference pixel [deg] |
| `CRPIX1`, `CRPIX2` | Reference pixel (by default the frame centre) |
| `CD1_1`, `CD1_2`, `CD2_1`, `CD2_2` | Scale, rotation and reflection matrix [deg/pix] |
| `CUNIT1`, `CUNIT2` | `deg` (`pix` for relative astrometry) |
| `CRDER1`, `CRDER2` | RMS of the residuals [deg] |
| `PV1_0…PV1_11`, `PV2_0…PV2_11` | TPV coefficients (`-p TPV` only) |
| `RADESYS`, `JEPOCH` | Reference system and epoch |
| `ASTRMETH` | How the match was found: `WCS`, `POINTING`, `BLIND`, `CROSS` or `SEQUENCE` |

A `COMMENT` block between `=== Astrometric Solution by Munipack ===` and
`=== End of Astrometric Solution by Munipack ===` records:

- the reference and the projection;
- the match sequence length and the number of objects used;
- the fitting method and the match type;
- the RMS, and the scale and rotation angle with their errors;
- the reflection and the projection centre.

With `--report` the extension **`ASTROREPORT`** lists the stars of the fit:

| Column | Description |
|---|---|
| RA, Dec (named as in the catalogue, e.g. `RAJ2000`, `DEJ2000`) | Catalogue coordinates [deg] |
| `X`, `Y` | Frame coordinates [pix] |
| `RA_OFF`, `DEC_OFF` | Residuals [deg] |
| `MATCH` | Position of the star in the match sequence; 0 for stars added by the final cross-match |

## Examples

Calibrate a frame from scratch against a cone around the target:

```bash
munipack cone -c UCAC5 -r 0.2 -o cone.fits -- 110.25 71.34
munipack find frame.fits
munipack aphot frame.fits
munipack astrometry -c cone.fits frame.fits
```

The pointing and the scale are known but the roll is not (common for
space-borne or alt-az instruments). This example uses a 0.756″/px detector
(4761.9 pix/deg) and a pointing known to a few arcminutes:

```bash
munipack astrometry -c cone.fits --rcen 180.02 --dcen 0.51 --scale 4761.9 \
    --wcs-offset 3arcmin frame.fits
```

The frames already have an approximate header WCS whose rotation may be
wrong. This skips the pass within 3°, writes the TPV distortion and the
residual table, and uses four threads (in the container, pass
`-e OMP_NUM_THREADS=4` to `docker run` instead):

```bash
OMP_NUM_THREADS=4 munipack astrometry --roll-unknown -p TPV --report \
    -c cone.fits -t calibrated/ frame_*.fits
```

Write a known calibration by hand (no catalogue), or remove one:

```bash
munipack astrometry -m manual --rcen 30 --dcen 60 --xcen 300 --ycen 300 \
    --scale 3600 --angle 20 frame.fits
munipack astrometry --remove frame.fits
```

```{note}
In `manual` mode, give `--xcen` and `--ycen` explicitly. If they are left
out, this build writes `CRPIX1 = CRPIX2 = 0` instead of the frame centre.
```

```{note}
Relative astrometry (`-R`) failed with "No common stars found" in tests of
this build on simulated frames that `-r` and `-c` solve. Use `-r` (a
calibrated reference frame) or a catalogue instead.
```
