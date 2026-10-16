# find

Detects stars on calibrated frames and stores them in a `FIND` table
extension of each frame. The algorithm is based on DAOFIND (P. B. Stetson,
DAOPHOT II). The frame is convolved with a lowered Gaussian of the given FWHM,
and local maxima of the convolution above a threshold become candidates. Two
shape statistics, *sharpness* and *roundness*, separate stars from hot pixels,
cosmic-ray hits and bad columns. A candidate in the wing of a brighter star,
or near a saturated area, must also stand above the light of that wing. The
sky level is a mean clipped by the noise model of the frame (`GAIN`,
`READNS` or `RDNOISE`, `BIASLEV`), and the threshold height is computed
from the same model, so it keeps the same false-alarm rate on a sky of a
few photons per pixel as on a bright ground-based sky. `aphot`, `psfphot` and `astrometry` read the `FIND` table.

## Synopsis

```bash
munipack find [options] file(s)
```

In the container:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
    find -f 3.5 -th 5 frame.fits
```

## Options

### Detection

| Option | Argument | Default | Description |
|---|---|---|---|
| `-f`, `--fwhm` | pixels | 3 | FWHM of the detection filter. Use the FWHM of the stars: a filter much narrower than the PSF finds more noise peaks and splits stars. |
| `-th`, `--threshold` | *t* (sigmas) | 7 | Detection threshold, given as the Normal deviate of a false-alarm probability (see {ref}`find-threshold`). |
| `--false-rate` | *N* | — | Sets the threshold from the expected number *N* of noise detections per frame. Overrides `-th`. Must be > 0. |
| `--lothresh` | sigmas | 5 | Lower data limit: pixels below sky − *u*σ are not used in the detection. |
| `--saturate` | ADU | header `SATURATE`; otherwise 2<sup>BITPIX</sup> − 1 for integer frames, 0.99 × the maximum pixel for floating-point frames | Saturation level. Pixels at or above it are excluded, and saturated areas are not reported as stars. |
| `--read-noise` | ADU | frame header (`READNS` in ADU, or `RDNOISE` in e⁻ divided by `GAIN`) | Read noise. In `find` it only enters the noise model used to clip the sky; the threshold uses the observed scatter of the sky. |

### Shape tests

| Option | Argument | Default | Description |
|---|---|---|---|
| `--shape-sigma` | *k* (sigmas) | 3 | Accepts sharpness and roundness within *k* standard deviations of the value expected for a star. With `0` only the fixed windows below are used. |
| `--shrplo` | value | 0.2 | Lower sharpness limit. Giving it (or `--shrphi`) makes the sharpness window fixed. |
| `--shrphi` | value | 1.0 | Upper sharpness limit (fixed window). |
| `--rndlo` | value | −1.0 | Lower roundness limit. Giving it (or `--rndhi`) makes the roundness window fixed. |
| `--rndhi` | value | 1.0 | Upper roundness limit (fixed window). |

### Output files and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Writes the output frames to this directory and leaves the input frames untouched. |
| `-O` | — | off | Turns on advanced output file names (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against input names; its groups are `\1`, `\2`, … (for `-O`). |
| `--mask` | string | — | Output name built from the groups of `--pattern`, e.g. `'\1_f.\2'` (for `-O`). |
| `--format` | printf format | — | Numbered output names, e.g. `'out_%03d.fits'` (for `-O`). |
| `--verbose` | — | off | Prints the sky, the number of stars and the threshold statistics of each frame. |
| `-h`, `--help` | — | — | Shows the help. |

Input files can also be given as a list: `@list.txt`, or `-` to read names
from standard input.

(find-threshold)=
## Threshold

The convolution *h* of the sky has a scatter σ<sub>h</sub>. The threshold *t*
(`-th`) is a **false-alarm probability**: a pixel of the convolved sky exceeds
the threshold height with probability *Q*(*t*) = erfc(*t*/√2)/2, the same as
a Normal deviate above *t*. For Normal noise the height is simply
*t*σ<sub>h</sub>. A faint sky of a few photo-electrons per pixel is Poisson,
so its convolution is skewed with a heavier upper tail, and *t*σ<sub>h</sub>
would let up to 5 times more noise peaks through. When the frame has `GAIN`,
`find` therefore computes the height from the noise model: Poisson photons of
λ = (sky − `BIASLEV`)·`GAIN` electrons plus Normal noise for the rest of the
variance, inverted by a saddle-point approximation. Some examples:

- On a sky of 0.5 e⁻/pixel, `-th 4` at FWHM 2.7 needs a height of
  4.6σ<sub>h</sub>.
- On ground-based frames (hundreds of e⁻ per pixel) the height stays within
  0.03σ<sub>h</sub> of *t*.
- Without `GAIN` the noise is taken as Normal (*t*σ<sub>h</sub>).

The convolution is tested in every pixel, so even ideal noise produces false
detections. On a 2048² frame (4.2·10⁶ pixels) you can expect:

| Setting | Noise detections per 2048² frame |
|---|---|
| `-th 4` | ≈ 130 |
| `--false-rate 10` (*t* ≈ 4.58) | 7–13 |
| `-th 5` | ≈ 1 |
| `--false-rate 1` (*t* ≈ 5.04) | 0.2–1.8 |
| `-th 7` (default) | ≈ 5·10⁻⁶ |

With `GAIN` in the header these counts hold for skies from 0.1 to 200 e⁻ per pixel (the value for `-th 7` is the expectation, not a measurement). The expected count
is the lower of *A Q*(*t*) (*A* is the number of pixels) and the expected
Euler characteristic of a Gaussian random field with the filter's
roughness. `--false-rate N` solves for the *t* at which this count equals
*N*, so *t* depends on the frame size: `--false-rate 1` gives *t* ≈ 5.04 on
2048² and a higher *t* on larger frames.

```{note}
The shape tests do not reject noise peaks: near the threshold they look like
faint stars. The false detections are set by the threshold alone, so pick
`-th` or `--false-rate` for the depth and contamination you need.
```

A stacked frame that still has the `GAIN` of a single exposure in its header
gets a threshold that is too high. Correct `GAIN` before running `find`.

## Shape tests

- **Sharpness** compares the central pixel with the mean of the other pixels
  in the kernel, relative to the height.
- **Roundness** is 2(*h*<sub>x</sub> − *h*<sub>y</sub>)/(*h*<sub>x</sub> +
  *h*<sub>y</sub>) of the heights of the marginal fits.

DAOPHOT's fixed windows are 0.2…1.0 for sharpness and −1…1 for roundness.
Photon noise moves faint stars out of these windows. With a low gain, 40% of
the stars at 4–5σ fall outside 0.2…1.0.

With `--shape-sigma k` (default 3) a candidate passes if it is inside the
fixed window, **or** within *k* predicted standard deviations of the value
expected for a star. The standard deviation comes from the first-order noise
of the pixels involved: the sky plus the photon noise of a point source of the
measured height. The expected value is the median over bright candidates
(≥ 10σ, at least 20 of them), or a Gaussian model if there are too few.
Sharpness is always limited to −1…2. The fixed windows are used alone:

- with `--shape-sigma 0`;
- when `--shrplo`/`--shrphi` (sharpness) or `--rndlo`/`--rndhi` (roundness)
  are given;
- for candidates inside the wing of a brighter star.

**Wing test.** A candidate within 8 FWHM (at least 20 px) of a brighter
candidate, or close to a saturated area, is accepted only if its convolution
exceeds the convolved wing at the same distance by the threshold times the
noise of the wing. This removes the secondary detections in haloes and around
saturated cores, and the twin detections that a too-narrow filter finds
inside a star's own profile.

## Output

`find` adds (or replaces) a binary table extension `EXTNAME = 'FIND'` in
each frame. The primary array and its header are not changed. Without `-t`
or `-O` the input file is updated in place.

Columns (sorted by `PEAKRATIO`, brightest first):

| Column | Description |
|---|---|
| `X`, `Y` | Position in pixels. |
| `PEAKRATIO` | Height of the convolved peak divided by the threshold height. |
| `SHARP` | Sharpness. |
| `ROUND` | Roundness. |

Header keywords of the `FIND` extension:

| Keyword | Description |
|---|---|
| `FWHM` | FWHM of the filter [pix]. |
| `THRESH` | *t*: the Normal deviate of the false-alarm probability (computed when `--false-rate` is used). |
| `HTHRESH` | The threshold height in σ<sub>h</sub> under the noise model of the sky. |
| `SKEWFILT`, `KURTFILT` | Skewness and excess kurtosis of the convolved sky. |
| `FALSERAT` | Expected number of noise peaks above the threshold per frame (the requested value with `--false-rate`). |
| `SHAPESIG` | *k* of the shape tests (0: fixed windows only). |
| `SHARPREF`, `ROUNDREF` | Expected sharpness and roundness of a star. |
| `SHRPLO`, `SHRPHI`, `RNDLO`, `RNDHI` | The fixed shape windows. |
| `LOWBAD`, `HIGHBAD` | Range of good data [ADU]. |
| `SATURATE`, `READNS` | Saturation level used, and the read noise of `--read-noise` or of the header `READNS` (0 without either) [ADU]. |
| `ECCENTR`, `INCL` | Mean eccentricity and inclination of the detected stars. |

Comment cards also record the sky level and its scatter. `--verbose` prints
the sky level, the number of stars, the threshold and the shape references
of each frame.

## Examples

Detect stars with the default threshold (only fairly bright stars, almost no
noise detections):

```bash
munipack find -f 3.5 frame.fits
```

A deeper search that allows about one noise detection per frame, printing the
threshold statistics. Run it on copies in `found/` and keep the originals:

```bash
mkdir -p found
munipack find -f 4.9 --false-rate 1 --verbose -t found frame_*.fits
```

Faint stars on a low-count frame (for example a space-UV exposure with a gain
of about 0.1 e⁻/ADU), with a fixed threshold and DAOPHOT's old fixed shape
windows for comparison:

```bash
munipack find -f 4.9 -th 4 frame.fits                    # noise-aware shape tests
munipack find -f 4.9 -th 4 --shape-sigma 0 frame.fits    # fixed windows only
```

Provide the saturation and read noise when the header lacks them, then go on
with aperture photometry:

```bash
munipack find -f 3 -th 5 --saturate 60000 --read-noise 3.7 frame.fits
munipack aphot frame.fits
```

```{note}
Frames made by `artificial` have a small watermark in the top-right corner,
which `find` detects as a few "stars". To avoid it, simulate with
`--no-watermark`, or drop detections inside the region given by `WMARKX1`,
`WMARKX2`, `WMARKY1`, `WMARKY2`.
```
