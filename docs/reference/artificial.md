# artificial — simulated frames

`artificial` makes artificial images of the sky that are photometrically
realistic. They are meant for testing reduction pipelines and for studying
observing effects when the truth is known. A frame is built in electrons from
a star field (a catalogue or a random field), optional elliptical galaxies, a
variable star and the sky background, all seen through a model telescope
(area, quantum efficiency, exposure) and a model PSF. Optionally the
atmosphere adds extinction and seeing that grow with airmass. Then come a
flat-field model, the dark current, the Poisson and read-out noise, cosmic
rays, the gain, the bias, and saturation. Every frame carries a WCS, the PSF it
was rendered with, and a table of the true star positions and fluxes. This
makes it suitable for any optical or UV imager, from a ground-based CCD to a
low-count space-UV camera.

## Synopsis

```bash
munipack artificial [options]
```

There is no input file. Frames are written to the current directory as
`artificial.fits` (one frame) or `artificial_001.fits`, `artificial_002.fits`, …
(with `--count`). Existing files of the same name are replaced.

```{note}
In the Docker image:
`docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes artificial ...`
```

## Options

### PSF

| Option | Argument | Default | Description |
|---|---|---|---|
| `--psf` | `SEEING`, `MOFFAT`, `GAUSS` | `SEEING` | Analytic PSF shared by all stars. `SEEING` is the diffraction pattern of the aperture (Airy disk) convolved with a Gaussian seeing spread: the most natural profile, and the slowest. `MOFFAT` is (1 + r²/α²)<sup>−β</sup>. `GAUSS` is exp(−r²/2σ²): fast, and suited to space-based imaging. Case-insensitive; other names are an error. |
| `--psf-file` | FITS file | — | Render stars with a tabulated PSF, constant or varying over the field. Overrides `--psf`. See {ref}`psf-files` below. |
| `--hwhm` | real [pix] | 1 | PSF width. **Despite its name this is not the HWHM.** It is the Gaussian σ for `GAUSS` and the σ of the seeing spread for `SEEING` (HWHM = 1.1774 σ). For `MOFFAT` it is α (HWHM = α √(2<sup>1/β</sup> − 1)). |
| `--beta` | real | 2 | Moffat exponent β (β > 0, β ≠ 1). |
| `--airy` | real [pix] | from `--diameter`/`--area` | Radius of the Airy disk, used by `SEEING`. By default it is 0.07/R px, where R [m] is the aperture radius (0.125 px for the default 1.12 m). |
| `--spread` | `AUTO`, `FFT`, `RANDOM` | `AUTO` | How `SEEING` convolves the diffraction pattern with the seeing. `FFT` is fast and smooth. `RANDOM` adds randomly shifted copies of the pattern: slower, slightly asymmetric, like a short exposure. `AUTO` uses `RANDOM` for `--exptime` < 3 s and `FFT` otherwise. |
| `--eccentricity` | 0 ≤ e < 1 | 0 | Elliptical PSF: the contour ellipse has semi-major axis a (the width above) and b = a √(1 − e²). Use it to model tracking errors. |
| `--inclination` | real [deg] | 0 | Angle of the major semi-axis from the x axis, counter-clockwise. |

### Telescope, detector and noise

| Option | Argument | Default | Description |
|---|---|---|---|
| `--exptime` | real [s] | 1 | Exposure time. |
| `--diameter` | real [m] | 1.12 | Aperture diameter. It sets the collecting area (unless `--area` is given) and the Airy radius (unless `--airy` is given). |
| `--area` | real [m²] | 0.985 (π·0.56²) | Collecting area. Without `--diameter` it also fixes the aperture radius used for the Airy disk. |
| `--qeff` | 0 … 1 | 1 | Quantum efficiency of the whole instrument (optics, filter and detector). |
| `--dark` | real [e⁻/s/pix] | 0 | Dark current. It adds dark·t electrons uniformly to the mean, with Poisson noise. The flat, the sky gradient, airmass, area and QE do not affect it. Must be ≥ 0. |
| `--readnoise` | real [e⁻ rms] | 0 | Gaussian read-out noise, added after the Poisson noise and before the gain. It is independent of the exposure time. Must be ≥ 0. |
| `--gain` | real [e⁻/ADU] | 1 | Conversion gain. Frames are written in ADU. |
| `--bias` | real [ADU] | 0 | Bias level (pedestal), added after the gain. Use a few σ<sub>R</sub>/gain with `--readnoise` and integer `-B`, so that negative values are not clipped. |
| `--saturate` | real [ADU] | the largest value of `-B` | Saturation, applied to the total value including the bias. By default it is 2<sup>B</sup> − 1 for unsigned integer frames (255 for `-B 8`, 65535 for `-B 16`, 2³¹ − 2 for `-B 32`) and 10³⁸ for floats. A larger value is lowered to that limit. |
| `--disable-noise` | — | off | Write the noiseless mean frame: no Poisson noise (stars, sky, dark) and no read noise. Cosmic rays and the watermark are still added. |
| `--gnoise` | real [ph/s/m²/arcsec²] | 0 | **Deprecated**; use `--dark` and `--readnoise`. Despite its name this is not Gaussian noise. It is an extra uniform background rate merged with the sky rate in quadrature, √(sky² + gnoise²). It is kept so that old simulations can be reproduced. |
| `--flat-model-file` | FITS file | — (unity) | Per-pixel sensitivity (flat field) that multiplies the light (sky, stars, galaxies) before the noise. Its size must match the frame. |
| `-B`, `--bitpix` | `-32`, `32`, `16`, `8` | `-32` | Output data type. Integer frames are unsigned (`-B 16` and `-B 32` by BZERO), rounded to integers and clipped at zero. |

### Cosmic rays

| Option | Argument | Default | Description |
|---|---|---|---|
| `--cosmics-freq` | real [1/s] | 0 | Mean rate of cosmic-ray events. Each frame gets a Poisson number of events with mean freq·t. |
| `--cosmics-cmax` | real [e⁻] | 1000 | Maximum intensity along a track. |
| `--cosmics-len` | real [pix] | 10 | Maximum track length. |

### Sky background

| Option | Argument | Default | Description |
|---|---|---|---|
| `--sky-mag` | real [mag/arcsec²] | 21 | Sky brightness at the frame centre (at zenith with `--atmosphere`), in the filter `-f`. Use `99` or fainter for no sky. |
| `--sky-grad-x`, `--sky-grad-y` | real [mag/arcsec²/pix] | 0 | Gradient of the sky brightness from the frame centre. Positive values make the sky fainter in that direction. Typical values are 0.0001–0.001 (moonlight, light pollution). Only the sky follows the gradient; the dark and `--gnoise` do not. |

### Stars and galaxies

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | FITS table | — | Catalogue of star positions and magnitudes, for example the output of `munipack cone`. Without it a random field is generated. |
| `--maglim` | integer [mag] | 13 | Random field only: the faintest magnitude. Each integer magnitude bin m = 0 … maglim holds Σ(m) = ⌊10⁻⁴ e<sup>1.1m</sup>⌋ stars per square degree (a rough fit to Gaia DR1) times the frame area. Down to 13 mag that is 239 stars/deg² in total: about 8 stars on the default 666² frame and 77 on 2048² at 1″/px. Large values on large fields use much memory. |
| `--col-ra` | column name | from `--vocat` | Right Ascension column [deg]. |
| `--col-dec` | column name | from `--vocat` | Declination column [deg]. |
| `--col-pm-ra`, `--col-pm-dec` | column name | from `--vocat` | Proper-motion columns. They must exist in the table, but stars are not moved: positions are used as tabulated. |
| `--col-epoch-ra`, `--col-epoch-dec` | column name | from `--vocat` | Reference-epoch columns. Accepted, but not used. |
| `--col-mag` | column name | `<filter>mag` with `-f` (e.g. `Vmag`), else from `--vocat` | Magnitude column, in the filter `-f` of `--photsys`. |
| `--vocat` | XML file | `VOcat_conf.xml` of the installation | Catalogue configuration that gives the default column names of known catalogues. |
| `--galcat` | FITS table | — | Catalogue of elliptical galaxies with Sérsic profiles (extension `ELLIPTICAL`: `RAJ2000`, `DEJ2000`, `MUEFF` [mag/arcsec²], `REFF` [deg], `PA` [deg], `E`, `SERSIC`). A template is `src/doc/artgalaxy.lst`; build the FITS file with `munipack fits --restore artgalaxy.lst`. |

### Photometric system

| Option | Argument | Default | Description |
|---|---|---|---|
| `-f`, `--filter` | name | `V` | Filter of the simulated frame. It selects the photon rate of the magnitudes and of the sky, and is written to `FILTER`. |
| `--photsys` | name | `Johnson` | Photometric system that defines the filter (e.g. `Johnson`, `Gunn` for SDSS ugriz). |
| `--phsystab` | FITS file | `photosystems.fits` of the installation | Table of photometric systems. |

### Geometry (WCS)

| Option | Argument | Default | Description |
|---|---|---|---|
| `--width`, `--height` | integer [pix] | 666, 666 | Frame size. Give both or neither. |
| `--rcen`, `--dcen` | real [deg] | 0, 0 | RA, Dec of the projection centre, placed at pixel CRPIX = NAXIS/2 (integer division). Give both or neither. |
| `--scale` | real [pix/deg] | 3600 (1″/px) | Plate scale, **in pixels per degree**. |
| `--fov` | real [deg] | — | Field of view along the frame **height**. If given, it overrides `--scale` (scale = height/fov). |
| `--angle` | real [deg] | 0 | Rotation of the frame around its centre, clockwise positive. |
| `-p`, `--projection` | WCSLIB code | `TAN` | Celestial projection. |

### Variable star (light curve)

A variable star is added only when one of `--lc-model`, `--lc-table` or
`--lc-fourier` is given. Give only one of them: if several are given,
`--lc-fourier` takes precedence over `--lc-table`, and `--lc-table` over
`--lc-model`. The magnitude is evaluated at mid-exposure (the `JD` keyword).

| Option | Argument | Default | Description |
|---|---|---|---|
| `--lc-model` | `wave`, `Cepheid`, `Algol` | — | Built-in curve (case-insensitive): a cosine wave, δ Cep, or Algol (β Per). Other names are an error. |
| `--lc-table` | FITS table | — | Time-series table (phase, relative flux, flux error), sorted by phase and interpolated by a smoothing spline. |
| `--lc-fourier` | FITS table | — | Fourier coefficients a<sub>n</sub>, b<sub>n</sub> of the relative flux over one period (two columns; row n holds a<sub>n−1</sub>, b<sub>n−1</sub>). The template `src/doc/fourier.lst` becomes FITS with `munipack fits --restore fourier.lst`. |
| `--lc-mag` | real [mag] | 15 | Magnitude zero point m₀ of the curve: m = m₀ + 2.5 log(A·F). |
| `--lc-amp` | real | 1 | Amplitude A. |
| `--lc-jd0` | real [JD] | 2451544.5 | Reference epoch; use the full JD, not MJD. |
| `--lc-per` | real [day] | 1 | Period (> 0). |
| `--lc-ra`, `--lc-dec` | real [deg] | `--rcen`, `--dcen` | Position of the variable. Give both or neither. |

### Atmosphere and time

| Option | Argument | Default | Description |
|---|---|---|---|
| `--atmosphere` | — | off | Model the atmosphere. The airmass X comes from the time, the site and the field centre. Stars are attenuated by e<sup>−k(X−1)</sup>. The sky becomes b·X·e<sup>−k(X−1)</sup>, where b is the zenith sky. The PSF width scales as X<sup>0.6</sup>. A frame whose field is below the horizon is skipped with a warning. |
| `--extk` | real | 0.1 | Extinction coefficient k (monochromatic); 0 turns off extinction. |
| `--long`, `--lat` | real [deg] | 0, 0 | Site longitude (+east) and latitude (+north). Give both or neither. |
| `--date` | `YYYY-MM-DD` | 2000-01-01 | Date of the first exposure start (UT). |
| `--time` | `HH:MM:SS` | 00:00:00 | Time of the first exposure start (UT). |
| `--count` | integer | 1 | Number of frames in the series. |
| `--timestep` | real [s] | 60 | Interval between the starts of consecutive frames. |

### Reproducibility, parallel runs and output

| Option | Argument | Default | Description |
|---|---|---|---|
| `--seed` | integer | random | Seed of the random generator. With it, every frame is reseeded from (seed, frame index), and the random star field from the seed. Frames are reproducible and do not depend on `--mproc`. Recorded in `ARTSEED`. |
| `-mproc`, `--mproc` | integer > 0 | 3 natively (`MUNIPACK_CAF_MPROC`), 1 in the Docker image | Number of parallel processes for a series (coarray build, run by `mpirun`), at most the number of logical CPUs. |
| `--Xproc` | string | `MUNIPACK_CAF_XPROC` | Extra arguments for the parallel launcher (e.g. `--oversubscribe`). |
| `--mask` | file name | `artificial.fits` / `artificial_???.fits` | Output file name. A run of `?` is replaced by the zero-padded frame index (`frame_??.fits` → `frame_01.fits`, …). Without `?`, every frame of a series overwrites the same file. |
| `--no-watermark` | — | off | Leave out the watermark. By default a 34×30 px box, 10 px from the top-right corner, is drawn 5√B electrons above the sky (B the mean sky in electrons per pixel, so it is invisible on frames with no sky), with a pattern cut out. Frames smaller than 44×40 px get no mark. `find` detects it as a few stars. |
| `--fits-keys-file` | text file | — | Extra header cards, one per line in the CFITSIO template syntax (`KEY = value / comment`). They override the cards written by `artificial`. |
| `--verbose` | — | off | Print the parameters and the progress. |
| `-h`, `--help` | — | — | Show the help and exit. |

## Detector noise model

The mean number of electrons in pixel (x, y) is

```{math}
\lambda(x,y) = F(x,y)\,\bigl[B(x,y) + S(x,y)\bigr] + D\,t
```

where S is stars and galaxies, F the flat model (`--flat-model-file`, unity by
default), D the dark current (`--dark`) and t the exposure time. The sky is

```{math}
B(x,y) = S_0\,10^{-0.4\,[g_x (x-x_0) + g_y (y-y_0)]} + G_0,
\qquad S_0 = b\,t\,A\,\eta\,s^2
```

Here b is the sky photon rate from `--sky-mag` (per s, m² and arcsec²), A the
area, η the QE, s the pixel size in arcsec, g<sub>x</sub>, g<sub>y</sub> the
gradients, (x₀, y₀) the frame centre (NAXIS/2) and G₀ the uniform part from
`--gnoise`. A star of photon rate f (from its magnitude in the filter `-f`)
contributes f·t·A·η electrons, spread by the PSF and integrated exactly over
the pixels at its sub-pixel position. With `--atmosphere`, η of stars and
galaxies is multiplied by e<sup>−k(X−1)</sup> and S₀ by X e<sup>−k(X−1)</sup>.
The recorded frame is

```{math}
I(x,y) = \frac{\mathrm{Po}\bigl(\lambda(x,y)\bigr) + N(0,\sigma_R)}{G} + P
```

Po is Poisson noise on all electrons, N(0, σ<sub>R</sub>) is the read noise
(`--readnoise`), G is the gain and P the bias. Cosmic rays and the watermark
are added in electrons before the division by the gain. The result is clipped
at zero (integer `-B` only) and at `--saturate`, and integer frames are
rounded. A background pixel thus has mean λ/G + P and variance
(λ + σ<sub>R</sub>²)/G² ADU² (before clipping). With `--disable-noise` the
frame is λ/G + P (plus cosmic rays and the watermark).

```{note}
Clipping at zero biases the mean of a low background with read noise in
integer frames. Add a `--bias` of several σ<sub>R</sub>/G.
```

(psf-files)=
## PSF files

`--psf-file` takes a FITS file whose primary header has `ZOOM`, the number of
PSF samples per pixel. Each image extension that has `RADCOO` (distance from the
field centre, deg) and `ANGCOO` (position angle, deg) is the PSF at that field
position. The PSF is interpolated linearly in angle around each ring and in
radius between rings. A file with just a 2-D primary image and `ZOOM` is a
single PSF for the whole field. The star centre is the sample at `CRPIX1/2` of
the extension, or else at NAXIS/2. The samples should sum to ZOOM² for the
whole flux; a warning is printed if the sum differs by more than 0.1%. The
field radius is computed from the frame centre and `--scale`, so the file must
match the frame's scale and centre (the same `--scale`, `--width` and
`--height` as the frame it was measured on). `psfphot --psf EMPIRICAL --psf-save`
writes this format; a file saved from a frame without a WCS has `RADCOO` in
pixels and needs `--scale 1`. The full specification is
`src/doc/dataform_artpsf.shtml` in the source tree.

## Outputs

Each output file has these HDUs:

| HDU | Content |
|---|---|
| primary | The image in ADU, BITPIX from `-B`. |
| `PSF` | The analytic PSF used, sampled with `ZOOM` (16) samples per pixel. Not written with `--psf-file`. |
| `CATALOGUE` | One row per star (the variable first): `X`, `Y` [pix, FITS 1-based] and `COUNT` [electron], the mean number of electrons of the star after the QE and the extinction but before the flat and the gain (divide by `GAIN` to compare with ADU). |

Primary header keywords:

| Keyword | Meaning |
|---|---|
| `DATE-OBS` | Exposure start (UT). |
| `JD` | Julian date at mid-exposure. |
| `EXPTIME`, `FILTER`, `AREA` | Exposure [s], filter, collecting area [m²]. |
| `GAIN`, `SATURATE` | [e⁻/ADU], [ADU]. |
| `DARKCUR`, `RDNOISE`, `BIASLEV` | `--dark` [e⁻/s/pix], `--readnoise` [e⁻], `--bias` [ADU]. |
| `ARTSEED` | `--seed` (only when given). |
| `PSF`, `HWHM`, `ECCENTR`, `INCL` | PSF name, width (σ, or α for Moffat, after seeing scaling), eccentricity, inclination. With `--psf-file` they still show the `--psf` settings and do not describe the file PSF (a `COMMENT` names the file). |
| `MOFFATB` / `AIRYDISK`, `SPREAD` | β for `MOFFAT`; Airy radius and spread method for `SEEING`. |
| `AIRMASS` | With `--atmosphere`. |
| `CTYPE1/2`, `CRPIX1/2`, `CRVAL1/2`, `CD1_1` … `CD2_2` | The WCS. |
| `LONGITUD`, `LATITUDE` | Site. |
| `WMARKX1`, `WMARKX2`, `WMARKY1`, `WMARKY2` | First and last pixel of the watermark (not written with `--no-watermark`). |

`COMMENT` cards record the sky level and dark level in electrons per pixel,
the gradient, the number of cosmic-ray events, the PSF or PSF file, and the
light-curve, galaxy and flat-model inputs.

## Examples

A random field on a 16-bit CCD-like detector, reproducible:

```bash
munipack artificial --width 1024 --height 1024 --maglim 16 \
    --psf GAUSS --hwhm 1.5 --exptime 60 --sky-mag 21 \
    --gain 2 --readnoise 5 --bias 100 --dark 0.01 -B 16 \
    --seed 42 --no-watermark --mask field.fits
```

A real field from a catalogue, with the truth for testing `find`, `aphot`,
`astrometry` and `phcal`:

```bash
munipack cone -r 0.2 -o cat.fits -- 180.0 30.0
munipack artificial -c cat.fits -f V --rcen 180 --dcen 30 --scale 3600 \
    --col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA --col-pm-dec pmDE --col-mag Vmag \
    --width 640 --height 640 --exptime 60 --psf GAUSS --hwhm 1.5 \
    --gain 2 --readnoise 5 --bias 100 --seed 2 --mask frame.fits
```

A three-hour time series of an eclipsing binary seen from a ground site, with
extinction and seeing growing with airmass (1.4 to 2.9), on four processes:

```bash
munipack artificial --lc-model Algol --lc-mag 12 --lc-amp 1 \
    --lc-per 2.867 --lc-jd0 2460000.42 --rcen 47.04 --dcen 40.96 \
    --atmosphere --extk 0.15 --long 16.6 --lat 49.2 \
    --date 2023-02-24 --time 20:00:00 --count 36 --timestep 300 \
    --psf MOFFAT --hwhm 2.0 --beta 2.5 --exptime 30 \
    --seed 1 --mproc 4 --mask 'algol_??.fits'
```

This writes `algol_01.fits` … `algol_36.fits`. Frames whose time puts the
field below the horizon would be skipped with a warning.

Re-render a field with the empirical PSF measured on a real 640 × 640 frame
at 1″/px. The geometry must be that of the frame, since the PSF file is
positioned by the field centre and the scale:

```bash
munipack find -f 4 real.fits
munipack psfphot --psf EMPIRICAL --psf-save epsf.fits real.fits
munipack artificial -c cat.fits --rcen 180 --dcen 30 --width 640 --height 640 \
    --scale 3600 --psf-file epsf.fits --exptime 30 --seed 8 --mask sim.fits
```
