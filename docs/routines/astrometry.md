# astrometry

`astrometry` finds the WCS of a frame, so that every pixel gets a sky position.

It pairs the stars detected on the frame with a reference
catalogue (or an already calibrated frame), fits a projection to those
pairs, and writes the result into the primary header. The image data are
not changed.

If the header already has a rough WCS, or you know where the telescope
pointed, the match uses that and is fast. Otherwise it searches blindly.
Each solution is checked against chance coincidences before it is saved.

## Usage

```bash
munipack astrometry [options] file(s)
munipack astrometry [options] 'input.fits output.fits'
```

Files are updated in place. To keep the original, give the input and output
name as one quoted argument, or write into a directory with `-t`. A list of
files can be read with `@list.txt`, or from standard input with `-`. If one
frame fails, the others are still processed and the exit status is 1.

**Before you start**, run {doc}`find` (required) and
{doc}`aphot` (strongly recommended) on the frame, and get a reference
catalogue, for example with {doc}`cone`. `astrometry` uses the photometry
to pick the brightest stars and to weight the fit. There is no default catalogue, so a
reference (`-c`, `-r` or `-R`) is always needed, except in `manual` mode.

## Modes

| Mode | Use it when |
|---|---|
| `match` (default) | Almost always. Uses the header WCS or your pointing if available, otherwise matches blindly. |
| `cross` | You know the calibration approximately (header WCS or `--scale`, `--angle`, `--reflex`) and just want the stars paired by position. |
| `sequence` | You pair the stars yourself with `--seq1` and `--seq2`. |
| `manual` | You want to write a known WCS by hand, with no catalogue. |

## Options

### Mode and reference

| Option | Default | Description |
|---|---|---|
| `-m`, `--mode` *mode* | `match` | `match`, `cross`, `sequence` or `manual`. |
| `-c`, `--cat` *file* | — | Reference catalogue (FITS table). FITSIO filters work, e.g. `'cone.fits[1][$f.mag$>12]'`. |
| `-r`, `--ref` *file* | — | Calibrated frame used as the reference. |
| `-R`, `--rel` *file* | — | Reference frame for relative astrometry (pixel coordinates, no projection). |
| `--vocat` *file* | built-in | VO catalogue configuration that sets the default column names. |
| `--col-ra` *label* | `RAJ2000` | Right Ascension column [deg]. |
| `--col-dec` *label* | `DEJ2000` | Declination column [deg]. |
| `--col-pm-ra` *label* | `pmRA` | Proper motion in RA column. |
| `--col-pm-dec` *label* | `pmDE` | Proper motion in Dec column. |
| `--col-epoch-ra` *label* | by catalogue | Epoch of RA column. |
| `--col-epoch-dec` *label* | by catalogue | Epoch of Dec column. |
| `--col-mag` *label* | by catalogue | Magnitude column, used to pick the brightest stars. |

Use only one of `-c`, `-r` and `-R`. Proper motions are applied to the
`DATE-OBS` of the frame when the catalogue gives an epoch.

### Projection and pointing

| Option | Default | Description |
|---|---|---|
| `-p`, `--projection` *name* | `TAN` | `TAN` (gnomonic) or `TPV` (gnomonic with polynomial distortion). |
| `--xcen`, `--ycen` *pix* | frame centre | Reference pixel. Give both. |
| `--rcen`, `--dcen` *deg* | — | Sky position of the reference pixel; in `match` mode, the pointing. Give both. |
| `--scale` *pix/deg* | — | Plate scale in **pixels per degree** (1″/px = 3600). |
| `--angle` *deg* | — | Rotation, clockwise positive. |
| `--reflex` | off | The frame is mirrored. |
| `--radesys` *name* | `ICRS` | Reference system written to the header. |
| `--epoch` *year* | from catalogue | Reference epoch written to the header. |

### Guided match

These tune how far the header WCS or your pointing may be off.

| Option | Default | Description |
|---|---|---|
| `--blind` | off | Ignore the header WCS and the pointing (`--scale` is still used). |
| `--scale-tol` *r* | 0.05 | Allowed relative error of the known scale; 0 means any scale. |
| `--scale-free` | off | Match with any scale (same as `--scale-tol 0`). |
| `--wcs-offset` *value* | 10% of frame (≥ 100 px) | Allowed pointing error, in pixels or with a unit: `90arcsec`, `3arcmin`, `0.1deg`, `20%`. |
| `--wcs-rot` *deg* | 3 | Allowed rotation error; 180 means any rotation. |
| `--wcs-scale` *r* | 0.03 | Allowed relative scale error of the header WCS. |
| `--roll-unknown` | off | The rotation is unknown: try any rotation straight away. |
| `--roll-known` | off | The rotation is known within `--wcs-rot`: don't try other rotations. |
| `--reflex-unknown` | off | Try both mirror orientations. |

### Matching

| Option | Default | Description |
|---|---|---|
| `--sig` *pix* | 1 | Position error of the frame stars. |
| `--sigcat` *deg* | 1″ | Position error of the reference stars. |
| `--flxmin`, `--flxmax` *r* | 0.63, 1.58 | Allowed flux ratio of paired stars (±0.5 mag). Widen it if your filter differs from the catalogue band. |
| `--minmatch` *n* | 5 | Fewest stars in a match (at least 5). Use 7 or more on crowded fields. |
| `--maxmatch` *n* | 42 | Stars per chunk compared at once. |
| `--maxstars` *n* | 250 | How many of the brightest stars are used. |
| `--luckymatch` *n* | 3 | Accept a match once it is this many stars longer than `--minmatch`. |
| `--disable-lucky-match` | off | Stop at the first match found. |
| `--disable-flux-check` | off | Don't compare fluxes while matching. |
| `--disable-rms-check` | off | Skip the quality checks of the solution. Use with care. |
| `--seq1` *i,j,…* | — | Reference star indices (`sequence` mode). |
| `--seq2` *i,j,…* | — | Matching frame star indices (`sequence` mode). |

### Fit and output

| Option | Default | Description |
|---|---|---|
| `--fit` *method* | `robust` | `robust` or `squares`. |
| `--units` *unit* | auto | Units of the RMS in the header comments: `deg`, `arcmin`, `arcsec`, `mas`, `uas` (`pix` is accepted but ignored). |
| `--disable-save` | off | Compute the solution but don't write it. |
| `--remove` | off | Remove the astrometric calibration from the header. |
| `--report` | off | Add an `ASTROREPORT` table of the stars used in the fit. |
| `--show-defaults` | — | Print the default parameters and exit. |
| `-t`, `--target-directory` *dir* | — | Write output files into this directory. |
| `-O` | off | Enable advanced output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` *regex* | `(.+)\.(.+)` | Pattern matched against input names. |
| `--mask` *string* | — | Output name from back-references `\1`, `\2`, … |
| `--format` *printf* | — | Numbered output names, e.g. `out_%02d.fits`. |
| `-mproc`, `--mproc` *n* | 1 | Number of parallel processes. |
| `--Xproc` *string* | — | Extra options for `mpirun`. |
| `--verbose` | off | Print details of the match and the fit. |
| `-h`, `--help` | — | Show help and exit. |

## Output

The WCS is written to the primary header (`CTYPE`, `CRVAL`, `CRPIX`, `CD`,
plus `PV` terms for TPV) and replaces any earlier solution. `CRDER1/2` hold
the RMS of the fit, and `ASTRMETH` records how the match was found
(`WCS`, `POINTING`, `BLIND`, `CROSS` or `SEQUENCE`). A `COMMENT` block
summarises the solution: reference, number of stars, RMS, scale and
rotation. With `--report`, the `ASTROREPORT` extension lists each star's
catalogue and frame coordinates and its residuals.

## Examples

Calibrate a frame from scratch:

```bash
munipack cone -c UCAC5 -r 0.2 -o cone.fits -- 110.25 71.34
munipack find frame.fits
munipack aphot frame.fits
munipack astrometry -c cone.fits frame.fits
```

You know the pointing and the pixel scale, but not the rotation:

```bash
munipack astrometry -c cone.fits --rcen 180.02 --dcen 0.51 --scale 4761.9 \
    --wcs-offset 3arcmin frame.fits
```

Write a known calibration by hand, or remove one:

```bash
munipack astrometry -m manual --rcen 30 --dcen 60 --xcen 300 --ycen 300 \
    --scale 3600 --angle 20 frame.fits
munipack astrometry --remove frame.fits
```

## Notes

**Catalogue size.** If the pointing is known, a cone larger than the frame
is fine, since only stars near the frame are used. For a blind match, make
the cone just a bit larger than the field of view (radius about 0.8 × the
frame side for a square frame).

**Magnitudes.** Choose a catalogue magnitude range that matches the stars
that are well exposed but not saturated on your frame, using
`cone --magmin/--magmax` or a FITSIO filter.

**Crowded fields and few stars.** At least five stars must be paired; seven
or more is safer, so raise `--minmatch` to 7 on crowded fields. With few
stars, a tight `--wcs-offset` helps, because a smaller search is less
likely to be fooled by chance alignments.

**Distortion.** `-p TPV` fits a polynomial distortion on top of the
gnomonic projection. It needs plenty of stars: with 25 matched stars or
fewer it falls back to `TAN` with a warning.

**Speed.** The blind match runs in parallel threads; set the number with
`OMP_NUM_THREADS`.

```{note}
In `manual` mode, always give `--xcen` and `--ycen`. Without them the
reference pixel is written as 0 rather than the frame centre.
```

Relative astrometry with `-R` is unreliable in this version; use a
calibrated reference frame (`-r`) or a catalogue instead.
