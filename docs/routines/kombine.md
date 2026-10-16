# kombine

`kombine` combines astrometrically calibrated frames into one image: a deep
stack of a series, a mosaic of several fields, or one frame resampled onto a
new grid.

Each frame is reprojected through its own WCS onto the output grid and the
frames are averaged pixel by pixel. By default the average is a mean clipped
by the noise model of the frames, so cosmic rays, hot pixels and satellite
trails are rejected. The stored image is N × the mean of N frames, which for
equal exposures is simply their sum. The background of each frame is
subtracted first.

## Usage

```bash
munipack kombine [options] file(s)
```

The inputs must be 2-D frames with a WCS, for example from {doc}`astrometry`.
Frames without one are skipped with an error message and the rest are still
combined. If you want to do photometry on the result, calibrate the inputs
with {doc}`phcorr` first.

## Options

Any grid parameter you leave out is taken from the first input frame. Options
listed as a pair work only when both are given.

### Output grid

| Option | Argument | Default | Description |
|---|---|---|---|
| `-p`, `--projection` | code | first frame's (usually `TAN`) | WCS projection code, e.g. `TAN` or `SIN`. |
| `--width`, `--height` | integer [pix] | first frame's size | Size of the output. Give both. |
| `--xcen`, `--ycen` | real [pix] | centre of the output | Pixel where the projection centre falls. Give both. |
| `--rcen`, `--dcen` | real [deg] | centre of the first frame | RA and Dec of the projection centre. Give both. |
| `--scale` | real [pix/deg] | first frame's scale | Output scale; 3600 means 1″ per pixel. |
| `--angle` | real [deg] | first frame's rotation | Position angle; `0` puts north up. |
| `--reflex` | `yes` \| `no` | first frame's | Mirror the output or not. |

### Moving targets

| Option | Argument | Default | Description |
|---|---|---|---|
| `--pm-ra`, `--pm-dec` | real [deg/day] | 0, 0 | Motion of the target, so that it stays sharp and the stars trail. Give both. |
| `--pm-jdref` | real [JD] | mid-exposure of the first frame | Time at which the target position is taken. |

### Resampling and averaging

| Option | Argument | Default | Description |
|---|---|---|---|
| `-i`, `--interpol` | `near`, `bilinear`, `bicubic`, `bi3conv` | `near` | Interpolation method, see {ref}`kombine-resampling`. |
| `--flux-conserve` | — | off | Conserve the total flux of objects instead of surface brightness. |
| `--disable-back` | — | off | Keep the background instead of subtracting it. |
| `--arith` | — | off | Plain exposure-weighted mean: no rejection, much less memory. |
| `--hampel` | — | off | The older Hampel robust mean, for reproducing old results. |

### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `kombine.fits` | Output file (overwritten if it exists). |
| `-B`, `--bitpix` | `-32`, `32`, `16`, `8` | `-32` | Data type of the output. |
| `--verbose` | — | off | Print a summary of the output and warnings. |
| `-h`, `--help` | — | — | Show the help and exit. |

```{warning}
Keep the default `-B -32`. With an integer BITPIX the error planes of the
output often do not fit, and `kombine` then stops without writing anything.
```

(kombine-resampling)=
## Resampling

- `near` (the default) is the fastest and keeps output pixels independent,
  but undersampled frames may look a little blocky.
- `bilinear` and `bicubic` give smoother images; `bi3conv` is meant for
  strongly undersampled frames. They can leave faint fringes on the
  background where frames overlap.

Interpolation mixes neighbouring pixels, so the output noise is correlated:
per-pixel errors underestimate the noise in an aperture. There is no
drizzle-style resampling, so avoid strong downsampling.

By default the output preserves surface brightness: the sky per pixel stays
the same, while star fluxes change with the pixel area. Halving the scale
(2″ pixels instead of 1″) makes star fluxes four times smaller. Use
`--flux-conserve` when you need total fluxes for photometry.

## Output

One FITS file with three HDUs:

- **primary**: the combined image (the sum for equal exposures), with the
  WCS of the output grid, the total `EXPTIME`, the earliest `DATE-OBS`, and
  `FILTER` and `OBJECT` from the first frame.
- **STDERR**: the standard error of the image, on the same scale; NaN where
  only one frame covers a pixel.
- **MAXDEV**: the largest deviation of each pixel's values from the mean, in
  units of σ; handy for spotting outliers. Not written with `--arith`.

## Examples

Stack a series onto the grid of the first frame:

```bash
munipack kombine -o stack.fits frame_*.fits
```

Build a 1000 × 1000 mosaic at 2″/pixel, north up, keeping total fluxes and
the sky for photometry:

```bash
munipack kombine --rcen 201.365 --dcen -43.019 --width 1000 --height 1000 \
    --scale 1800 --angle 0 --reflex no --flux-conserve --disable-back \
    -i bilinear -o mosaic.fits field_*.fits
```

Follow a comet moving 0.4°/day in RA and −0.01°/day in Dec:

```bash
munipack kombine --pm-ra 0.4 --pm-dec -0.01 -o comet.fits comet_*.fits
```

## Notes

- The default and `--hampel` keep every reprojected frame in memory (about
  400 MB for 100 frames of 1000 × 1000 pixels). For very large sets of frames
  already free of cosmic rays, use `--arith`. The two options cannot be
  combined.
- The clipped mean reads `GAIN`, `READNS` (or `RDNOISE`) and `BIASLEV` from
  each frame. Without `GAIN` it relies on the scatter of the values only.
- The RA rate of `--pm-ra` is the rate of the coordinate itself, not
  multiplied by cos δ.
- To find the scale of a frame, read its CD matrix
  (`munipack fits -K CD1_1,CD1_2 frame.fits`); the scale in pix/deg is
  1/√(CD1_1² + CD1_2²).
- To make a colour image, resample each band with `kombine` onto the same
  grid and pass the results to {doc}`colouring`.
