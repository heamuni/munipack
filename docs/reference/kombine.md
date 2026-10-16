# kombine, colouring — stacking, mosaics and colour frames

## kombine

`kombine` combines frames that already have an astrometric (WCS) calibration
into a single image: a deep stack, a mosaic, or a frame resampled to a new
grid. You can choose the output grid (size, centre, scale, rotation,
reflection, projection). Any parameter you leave out is taken from the first
input frame. Each input frame is reprojected onto the grid through its own
WCS. The value of every input pixel is interpolated at the output pixel
centres and divided by its relative exposure, w = EXPTIME / mean EXPTIME.
`kombine` then averages the values of each output pixel. By default it uses a
mean clipped by the noise model of the frames. This mean rejects cosmic rays,
hot pixels and satellite trails, and stays unbiased on a sky of a fraction of
an electron per pixel. The stored image is N × the mean, where N is the
number of frames. For equal exposures that is the sum of the frames. By
default each frame's background is subtracted before averaging.

### Synopsis

```bash
munipack kombine [options] file(s)
```

The input frames must be 2-D and carry a WCS, for example from
[`astrometry`](astrometry.md). A file that cannot be read, or has no
astrometric calibration, is skipped with an error message. The remaining
files are still combined. To use the result for photometry, calibrate the
inputs first (`phcorr`).

In a container:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
    kombine -o stack.fits frame_*.fits
```

### Options

#### Output grid (projection)

Options marked as a pair take effect only when both are given. If one of them
is missing, the pair is ignored.

| Option | Argument | Default | Description |
|---|---|---|---|
| `-p`, `--projection` | code | the first frame's (TAN for nearly all frames) | Spherical projection of the output, as a three-letter WCS code (case-insensitive), e.g. `TAN` (gnomonic) or `SIN`. |
| `--width`, `--height` | integer [pix] | the first frame's size | Size of the output frame. Give both. |
| `--xcen`, `--ycen` | real [pix] | the output centre, (width + 1)/2, (height + 1)/2 | Output pixel (CRPIX) that the projection centre falls on. Give both. |
| `--rcen`, `--dcen` | real [deg] | the sky position of the first frame's centre | Projection centre (CRVAL) in Right Ascension and Declination. Give both. |
| `--scale` | real [pix/deg] | the first frame's scale | Output scale in **pixels per degree**. For example, 3600 is 1″/pixel. |
| `--angle` | real [deg] | the first frame's rotation | Position angle of the output. `--angle 0` gives north up. |
| `--reflex` | `yes` \| `no` | the first frame's | Mirror the output (`yes`) or not (`no`). Any value containing `y` counts as yes. |

#### Moving targets

| Option | Argument | Default | Description |
|---|---|---|---|
| `--pm-ra`, `--pm-dec` | real [deg/day] | 0, 0 | Motion of the target in RA and Dec. Each frame is shifted by the motion between its mid-exposure and the reference time, so the target stays fixed on the output and the stars trail. The RA rate is dα/dt of the coordinate itself, not multiplied by cos δ. Give both. |
| `--pm-jdref` | real [JD] | the mid-exposure JD of the first frame | Reference time of the motion: the target is placed at its position at this time. |

#### Resampling and averaging

| Option | Argument | Default | Description |
|---|---|---|---|
| `-i`, `--interpol` | `near`, `bilinear`, `bicubic`, `bi3conv` | `near` | Interpolation method (case-insensitive): nearest neighbour, bilinear, bicubic, or bicubic convolution. Any other value is an error. See {ref}`kombine-resampling`. |
| `--flux-conserve` | — | off | Multiply each value by the area of the output pixel in input pixels (the local Jacobian). This conserves the total flux of objects. Without it, `kombine` conserves the surface brightness. |
| `--disable-back` | — | off | Keep the background. By default the robust mean of each frame is subtracted before averaging. |
| `--arith` | — | off | Use the exposure-weighted arithmetic mean instead of the clipped mean. Cannot be combined with `--hampel`. |
| `--hampel` | — | off | Use the robust mean of the Hampel M-estimator, the default before 0.6.2+hg1653. Cannot be combined with `--arith`. |

#### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `kombine.fits` | Output file. An existing file is overwritten. |
| `-B`, `--bitpix` | `-32`, `32`, `16`, `8` | `-32` | BITPIX of all output HDUs. With a positive BITPIX only IMAGE is clipped to 0 … 2<sup>BITPIX</sup> − 1. STDERR and MAXDEV are not clipped, and a value that does not fit stops `kombine` with a FITSIO "datatype conversion overflow" error and no output. With the default mean this happens on typical frames, so keep `-32`. |
| `--verbose` | — | off | Print the output size, object, filter, date and total exposure, plus warnings about missing keywords. |
| `-h`, `--help` | — | — | Show the help and exit. |

(kombine-resampling)=
### Resampling

The output grid is sampled at pixel centres, using the FITS convention.
Bilinear, bicubic and bicubic convolution fall back to bilinear or nearest
neighbour near the frame edges, where their stencil would leave the frame.
Output pixels that no frame covers are filled with the robust mean of the
covered pixels. That is about zero after background subtraction, or the mean sky with
`--disable-back`.

- `near` (the default) is the fastest. It keeps the output pixels independent,
  but undersampled frames look slightly blocky.
- `bilinear` and `bicubic` give smoother images. `bi3conv` is meant for
  critically undersampled frames. All three bi-methods can leave faint
  fringing patterns on the background, especially where frames with
  different backgrounds overlap.
- Interpolation averages neighbouring input pixels, so the output errors are
  correlated. The per-pixel variance drops to about 0.52× (bilinear) or
  0.70× (bicubic) at random sub-pixel shifts, and STDERR follows it. The
  noise of sums over many pixels hardly changes: 0.87× and 0.94× over
  5 × 5 pixels. Per-pixel errors therefore underestimate the noise of an
  aperture (Casertano et al. 2000, AJ 120, 2747).
- `kombine` point-samples the input, so strong downsampling aliases. It does
  not do area-weighted (drizzle) resampling.

**Surface brightness versus flux.** By default the output holds surface
brightness: the sky level per pixel is preserved, and the total flux of a
star scales with the output pixel area by (s<sub>in</sub>/s<sub>out</sub>)²,
where s is the scale in pixels per degree. For example, halving the scale
(`--scale 1800` from 3600 pix/deg, i.e. 2″ pixels from 1″ pixels) makes star
fluxes four times smaller. `--flux-conserve` keeps the total fluxes instead,
and the sky per pixel then scales by the same ratio. Both modes give the same
result when the output scale equals the input scale.

### Means

All modes average v = value / w, where w is the exposure relative to the mean
one, using only the frames that cover the pixel. In a mosaic the number of
those frames, m, is the coverage of the pixel.

| Mode | Option | Behaviour |
|---|---|---|
| Clipped mean | default | Keeps the values within *t* − 5σ … *t* + *q* + √(*q*² + (5σ)²), iterated from the median. *t* is the mean. σ is the larger of the observed scatter and the frame's noise model, √((*t* + background − BIASLEV)/GAIN + READNS²) [ADU]. *q* = 25/6 photo-electrons. Read from each frame: GAIN (1 when BUNIT = COUNTS), READNS [ADU] or RDNOISE [e⁻], BIASLEV (0 when absent). Without GAIN only the scatter is used. Unbiased down to fractions of an electron per pixel, and rejects cosmic rays. This is the same routine as the sky of `aphot` and the master of `dark`. |
| Arithmetic | `--arith` | Exposure-weighted mean of the rates. Exact, but cosmic rays and satellites pass straight through. Needs much less memory (see below). |
| Hampel | `--hampel` | The old robust mean. It is biased low on a sky of a few photons per pixel (−13% at 0.39 e⁻/px, −2% at 2 e⁻, −0.2% at 20 e⁻) and biased differently for stars and sky. Kept for reproducing old results. |

The default and `--hampel` keep every reprojected frame in memory, about
width × height × 4 bytes per frame. For example, 100 frames of 1000 × 1000
pixels take 400 MB. `--arith` keeps only running sums.

### Outputs

One FITS file (`kombine.fits` by default):

| HDU | Content |
|---|---|
| primary | IMAGE: N × the mean of every pixel (the sum for equal exposures). |
| `STDERR` | The standard error of IMAGE on the same scale (the sum): N √(c/Σw), with c = Σw(v − v̄)²/(m − 1). NaN where fewer than two values were used (one frame, or a pixel covered once). |
| `MAXDEV` | The default and `--hampel` only: the largest deviation of a pixel's values from the mean, in units of σ (±Inf where σ = 0). Use it to find outliers and anomalies. |

Header of the primary HDU:

- WCS of the output grid: CTYPE, CUNIT, CRVAL, CRPIX, CD, RADESYS.
- `EXPTIME`: the total exposure time of all frames.
- `DATE-OBS`: the start of the earliest frame.
- `FILTER` and `OBJECT`, taken from the first frame. A warning is printed if
  another frame differs.
- `LONGITUD` and `LATITUDE`, when the first frame has them.
- `COMMENT` lines with the input files, the mean used (clipped, arithmetic or
  robust Hampel), a note that STDERR is the error of IMAGE (the sum), and the
  interpolation method.

The noise-model keywords (GAIN, READNS/RDNOISE, BIASLEV) are not copied to
the output.

```{note}
STDERR changed in this fork (BUGS S20). It is now the error of IMAGE, the
sum, and is NaN for a single value. Earlier versions wrote the error of the
mean, which was N times too small, and `--arith` wrote Σimg² where the
scatter was not positive.
```

### Examples

Stack a series onto the grid of the first frame:

```bash
munipack kombine -o stack.fits frame_*.fits
```

Build a 1000 × 1000 mosaic, north up and not mirrored, at 2″/pixel
(1800 pix/deg). `--flux-conserve` keeps the total fluxes of stars, and
`--disable-back` keeps the sky (per output pixel it scales with the pixel
area, see above), for photometry:

```bash
munipack kombine --rcen 201.365 --dcen -43.019 --width 1000 --height 1000 \
    --scale 1800 --angle 0 --reflex no --flux-conserve --disable-back \
    -i bilinear -o mosaic.fits field_*.fits
```

To find the input scale, read the CD matrix of a frame
(`munipack fits -K CD1_1,CD1_2 frame.fits`). The scale in pix/deg is
1/√(CD1_1² + CD1_2²).

Track a comet moving at +0.4°/day in RA and −0.01°/day in Dec. The comet
stays sharp and the stars trail:

```bash
munipack kombine --pm-ra 0.4 --pm-dec -0.01 -o comet.fits comet_*.fits
```

Use the plain mean, for example on long-exposure frames already cleaned of
cosmic rays, to save memory on a large set:

```bash
munipack kombine --arith -o stack.fits frame_*.fits
```

## colouring

`colouring` builds a colour FITS file from frames taken through standard
photometric filters (typically Johnson B, V, R). It does this by applying a
linear transformation into a colour space such as CIE 1931 XYZ, so the image
approximates natural colours. The result is a 3-D cube with one plane per
band. If the input frames are photometrically calibrated, or you give their
ctph with `-q`, the output is in absolute intensities (eV/s/m²/arcsec²).
Otherwise the white balance comes from the weights (`-w`), a white spot or a
white star. Without any of these `colouring` warns that the colours are not
balanced ("weights are undefined"). All bands must be on the same pixel grid. Resample them with `kombine`
first, using the same `--rcen`, `--dcen`, `--width`, `--height`, `--scale`,
`--angle` and `--reflex`.

### Synopsis

```bash
munipack colouring [options] files
```

List the frames from short to long wavelength (B, V, R). That is the
opposite of the order of the XYZ components.

### Options

#### Colour space

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cspace-input` | name | guessed from FILTER | Colour space of the inputs, e.g. `'Johnson BVR'`. Three frames with FILTER `B`, `V`, `R` (exactly, no aliases) are taken as `'Johnson BVR'`. Frames `Z`, `Y`, `X` are taken as `'CIE 1931 XYZ'`. Other sets need `-c`, otherwise `colouring` stops with an error. |
| `--cspace-output` | name | `'CIE 1931 XYZ'` | Colour space of the output. `'Johnson BVR'` keeps the input bands, for further processing. |
| `--list` | — | — | Print the available colour-space transformations and exit. |
| `--ctable` | file | `ctable.dat` of the installation | Definition of the colour-transformation matrices. |
| `--phsystab` | file | `photosystems.fits` of the installation | Photometric systems table, which gives the effective wavelengths of the bands. |

#### Bands and white balance

Lists are separated by commas. Where the locale's decimal point is a comma,
use semicolons instead. A list must have one item per input file.

| Option | Argument | Default | Description |
|---|---|---|---|
| `-w`, `--weights` | `w1,w2,…` | — | Positive weights of the bands. |
| `-q`, `--ctphs` | `ctph1,ctph2,…` | from calibrated frames | Reciprocal quantum efficiency of each band (from photometric calibration). Gives intensity-calibrated output. |
| `-b`, `--backs` | `b1,b2,…` | robust mean of each frame | Background level of each band. |
| `--disable-back` | — | off | Do not estimate or subtract backgrounds. Useful for non-astronomical images. |
| `--white-spot` | `x,y` [pix] | — | Centre of a region that should come out white. The weights come from its robust mean in each band. |
| `--white-star` | `x,y` [pix] | — | Centre of a star taken as white (B − V = 0, an A0 star). The weights come from its aperture flux in each band. |
| `--white-radius` | real [pix] | 7 | Radius of the white spot or star aperture. |

#### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `colour.fits` | Output file. An existing file is overwritten. |
| `-B`, `--bitpix` | `-32`, `32`, `16`, `8` | `-32` | BITPIX of the output. |
| `--verbose` | — | off | Print the guessed colour space and progress. |
| `-h`, `--help` | — | — | Show the help and exit. |

### Outputs

A colour FITS file: a cube of NAXIS1 × NAXIS2 × the number of bands. The
header has:

- `CSPACE`: the output colour space. A `COMMENT` gives the input colour space.
- `COMMENT` lines with each file, its weight and its background.
- `OBJECT`, `DATE-OBS`, `EXPTIME` and `AREA`, taken from the V frame.
- The WCS of the V frame, when present, with the spectral axis `CTYPE3 = 'WAVE'`.
- `BUNIT = 'eV/s/m2/arcsec2'` for calibrated or ctph input (`pixel2` without astrometry).
- `CREATOR`.

### Examples

```bash
# list the available transformations
munipack colouring --list

# B, V, R frames on a common grid -> CIE 1931 XYZ colour frame
munipack colouring -o m27.fits m27_B.fits m27_V.fits m27_R.fits

# uncalibrated frames: white balance by a star of B - V = 0 at (512, 488)
munipack colouring --white-star 512,488 --white-radius 5 -o m27.fits \
    m27_B.fits m27_V.fits m27_R.fits
```
