# bias, dark, flat, phcorr — calibration frames and corrections

These four actions turn raw CCD or CMOS exposures into calibrated frames in
electrons:

| Action | Input | Output |
|---|---|---|
| `bias` | zero-exposure frames | master bias (per-pixel robust mean) |
| `dark` | dark frames, optionally a master bias | master dark (per-pixel mean clipped by the noise model) |
| `flat` | flat-field frames, optionally bias and dark | master flat, normalised to a mean of 1 |
| `phcorr` | science frames, any of the masters above | corrected frames in electrons, with a `STDERR` extension |

The masters are single files. `phcorr` writes one output per input frame.
The calibrated frames are then ready for `find`, `aphot`, `psfphot` and
`astrometry`.

In the container, mount the working directory on `/data`:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
    bias -o bias.fits bias_*.fits
```

(calibration-sequence)=
## Typical calibration sequence

```bash
# 1. master bias from many zero-exposure frames (more than 10 is better)
munipack bias -o bias.fits bias_*.fits

# 2. master dark, bias subtracted (so that it can be scaled by exposure time)
munipack dark -bias bias.fits -o dark120.fits dark120_*.fits

# 3. master flat for each filter, from bias- and dark-corrected flats
munipack flat -bias bias.fits -dark dark120.fits -o flat_V.fits flat_V_*.fits
munipack flat -bias bias.fits -dark dark120.fits -o flat_R.fits flat_R_*.fits

# 4. the science frames: gain, bias, scaled dark, flat; originals untouched
mkdir -p cal
munipack phcorr -bias bias.fits -dark dark120.fits -flat flat_V.fits \
    -t cal object_V_*.fits
```

`phcorr` computes, for every pixel,

*C* = *g* (*S* − *t D* − *B*) / *F*,

where *S* is the raw frame, *g* the gain, *B* the bias, *D* the dark, *t* the
dark time factor and *F* the flat-field (see [phcorr](#phcorr)). There are
two common strategies:

- **Darks of the science exposure time, no bias.** Use `dark` without
  `-bias` and `phcorr -dark` without `-bias`. The dark then contains the
  bias, and *t* = 1 because the exposure times are equal. This is the
  simplest choice when all frames have one exposure time.
- **Bias plus scaled darks.** Make the dark with `-bias`, then give the
  *same* bias to `flat` and `phcorr`. *t* is the ratio of the exposure
  times, so one long dark serves every exposure time.

```{note}
Do not scale a dark that still contains the bias. With a dark made without
`-bias`, an exposure-time ratio other than 1 scales the bias pedestal too.
Use such a dark only for frames of its own exposure time, or make the dark
with `-bias` and give the same bias to `flat` and `phcorr`.
```

## Behaviour common to bias, dark and flat

**Pixel rejection.** Input pixel values outside `threshold < value <
saturation` are left out of the average. Both limits are in the units of the
raw input frames (ADU) and apply before any correction.

- **Saturation (`-st`).** Defaults to the `SATURATE` header keyword. Without
  it, integer frames use 2<sup>BITPIX</sup> − 1 and float frames the type
  maximum, which turns the check off.
- **Threshold (`-th`).** Defaults to 1 for `bias` and `dark`, and to the
  machine epsilon for `flat`. It must be > 0, so frames with no pedestal and
  values near zero lose their low pixels.

**Masks (`-bitmask`).** A FITS image of the same size. Pixels with a value
above 0.5 are used and pixels ≤ 0.5 are skipped (`BITPIX = 8` is enough).
`bias` and `dark` set the masked pixels, and pixels with fewer than two valid
values, to the mean level of the master. `flat` leaves masked pixels out of
the fit and out of its normalisation, and sets them to 1 in the master.

**Header checks.** `bias` and `dark` compare the exposure time and chip
temperature of the inputs. A frame differing from the mean by more than
`MUNIPACK_EXPTIME_TOLERANCE` (default 10⁻⁶ s) or
`MUNIPACK_TEMPERATURE_TOLERANCE` (default 1 °C) gets a warning. `flat`
warns when the inputs have different or empty filters. If the keywords
have other names in your files, rename them with the `FITS_KEY_*`
environment variables. For example, `FITS_KEY_EXPTIME=EXPOSURE` reads the
exposure time from `EXPOSURE` (and writes it there in the master). Others
are `FITS_KEY_FILTER`, `FITS_KEY_TEMPERATURE`, `FITS_KEY_DATEOBS`,
`FITS_KEY_GAIN` and `FITS_KEY_SATURATE`. The read noise is always read from
`READNS` [ADU] or `RDNOISE` [e⁻]. In Docker, pass the variables with
`-e FITS_KEY_EXPTIME=EXPOSURE`.

**Output.** The default output type is `-B -32` (float), which is the one to
use for masters. The default output names are `bias.fits`, `dark.fits` and
`flat.fits`, written to the current directory. An existing file of that name
is replaced.

**Input lists.** Files can be given as arguments, as `@list.txt` (one name
per line), or as `-` to read names from standard input (in Docker, add `-i`
to `docker run`).

(bias)=
## bias

Averages zero-exposure (bias) frames into a master bias, which maps the
electronic offset of the readout. Each output pixel is the robust mean
(Hampel M-estimator) of that pixel over the input frames. Bias frames contain
only read noise, which is symmetric, so the robust mean is unbiased. It also
rejects outliers such as cosmic-ray hits. Bias exposures are short, so
average many of them: more than 10.

### Synopsis

```bash
munipack bias [options] file(s)
```

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `bias.fits` | Output file. |
| `-st`, `--saturate` | ADU | header `SATURATE` | Upper limit: pixels at or above it are not averaged. |
| `-th`, `--threshold` | ADU | 1 | Lower limit: pixels at or below it are not averaged. Must be > 0. |
| `-bitmask` | file | — | Mask frame: pixels ≤ 0.5 are skipped and set to the frame mean. |
| `-B`, `--bitpix` | −32, 32, 16, 8 | −32 | Numerical type of the output. |
| `--verbose` | — | off | Prints the mean level, exposure time and temperature. |
| `-h`, `--help` | — | — | Shows the help. |

### Outputs

- **Primary HDU:** the master bias, in the units of the input (ADU).
- **`STDERR` extension:** the standard error of the mean of each pixel.
- **Header:** `DATE-OBS` of the first input frame. `EXPTIME` and `TEMPERAT`:
  the averages of the inputs. `IMAGETYP` copied from the inputs.
  `DARKMEAN = 'HAMPEL'` (the averaging method). `COMMENT` lines with the mean
  level and the list of input files.

### Examples

```bash
# master bias from all bias frames in the directory
munipack bias -o bias.fits bias_*.fits

# frames listed in a file, with a bad-column mask
ls raw/bias*.fits > bias.lst
munipack bias -bitmask badcols.fits -o bias.fits @bias.lst

# 16-bit camera saturating at 60000 ADU
munipack bias -st 60000 -o bias.fits bias_*.fits
```

(dark)=
## dark

Averages dark frames into a master dark. The dark current is proportional to
the exposure time and depends steeply on the chip temperature. With `-bias`,
the master bias is subtracted from every frame before averaging, so the
master contains the dark current only. `phcorr` and `flat` can then scale it
to other exposure times.

By default each pixel is the **mean clipped by the noise model** of its
frames. The mean is taken over the values within *t* − 5σ … *t* + *q* +
√(*q*² + (5σ)²), iterated from the median. Here *t* is the mean, σ the larger
of the observed scatter and the model noise √((*t* − *b*)/`GAIN` + RN²), and
*q* = 25/6 photo-electrons.

- **RN** is the read noise, from `READNS` [ADU] or `RDNOISE` [e⁻].
- **b** is the zero-electron level: 0 with `-bias`, otherwise the header
  `BIASLEV`, otherwise 0.
- **Without `GAIN`**, σ is the observed scatter alone and *q* = 25/6 ADU.

This rejects cosmic rays and hot outliers. It does not clip the skewed Poisson
tail of a weak dark of a fraction of an electron per pixel. The robust mean
used before (`--hampel`) is biased low on such darks: about −18% at
0.13 e⁻/pixel. The clipped mean is unbiased from 0.05 to 20 e⁻/pixel.

### Synopsis

```bash
munipack dark [-bias file] [--arith | --hampel] [options] file(s)
```

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `-bias` | file | — | Master bias subtracted from every dark frame before averaging. Recorded as `DARKBIAS`. |
| `--arith` | — | off | Plain arithmetic mean of the pixel values, with no rejection. |
| `--hampel` | — | off | Robust mean (Hampel M-estimator) of earlier versions. Biased low on darks of a few electrons per pixel. Cannot be combined with `--arith`. |
| `-o`, `--output` | file | `dark.fits` | Output file. |
| `-st`, `--saturate` | ADU | header `SATURATE` | Upper limit of the raw values used. |
| `-th`, `--threshold` | ADU | 1 | Lower limit of the raw values used (before the bias is subtracted). Must be > 0. |
| `-bitmask` | file | — | Mask frame: pixels ≤ 0.5 are skipped and set to the frame mean. |
| `-B`, `--bitpix` | −32, 32, 16, 8 | −32 | Numerical type of the output. Keep −32 for a bias-subtracted dark, whose values can be near or below zero. |
| `--verbose` | — | off | Prints the noise model of the clipping (gain, read noise, zero-electron level), the mean level, exposure time and temperature. |
| `-h`, `--help` | — | — | Shows the help. |

### Outputs

- **Primary HDU:** the master dark (ADU), bias-subtracted when `-bias` is
  given. This also applies to a single input frame.
- **`STDERR` extension:** the error of the mean of each pixel,
  √(Σ(*x* − ⟨*x*⟩)²/(*n* − 1)/*n*), from the values kept. The error of the
  subtracted bias is not included. For a single input frame `STDERR` is zero.
- **Header:**
  - `EXPTIME`, `TEMPERAT`: the averages of the inputs.
  - `DATE-OBS`: from the first input frame.
  - `IMAGETYP`: copied from the inputs.
  - `DARKMEAN`: `CLIP`, `ARITH` or `HAMPEL`.
  - `DARKBIAS`: the name of the subtracted bias, with `-bias`.
  - `COMMENT` lines: the inputs, the bias and mask used, `-th` and `-st`.

  `phcorr` reads `EXPTIME` to scale the dark. If it reads `DARKBIAS` and
  is given the same bias, it propagates the bias error correctly.

### Examples

```bash
# darks of the science exposure time, bias left in (t = 1 in phcorr)
munipack dark -o dark120.fits dark120_*.fits

# bias-subtracted dark, to be scaled to any exposure time
munipack dark -bias bias.fits -o dark300.fits dark300_*.fits

# the exposure time is stored as EXPOSURE in these files
FITS_KEY_EXPTIME=EXPOSURE munipack dark -bias bias.fits -o dark.fits d_*.fits

# the same in the container
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data -e FITS_KEY_EXPTIME=EXPOSURE \
    munipack:0.6.2-hg1653-fixes dark -bias bias.fits -o dark.fits d_*.fits
```

(flat)=
## flat

Builds a master flat-field from flat frames: twilight sky, dome or night-sky
exposures. The flat maps the relative response of the whole light path:
pixel sensitivity, vignetting, dust shadows and filter non-uniformity. Each
input is first corrected as *g*(*F* − *t D* − *B*). Here *t* is the
exposure-time ratio of the flat to the dark, or `-xdark`, or 1 when no
exposure times are known. The frames are then scaled and averaged.

The standard method uses the Poisson statistics of the photons, where the
variance equals the mean in electrons. It solves pixel by pixel for the flat
that best fits all the scaled frames, iterating until the residuals fall
below the photon noise. Cosmic rays and other defects of single frames are
suppressed. The result is normalised to a mean of 1, and the filter is
copied from the inputs.

Guidelines:

- Use flats exposed to about half of the full well.
- Use more than 13 frames. The accuracy grows as √*N*.
- Memory use is about width × height × 8 bytes per frame: 100 frames of
  1000 × 1000 take 800 MB.
- Give each filter its own flat, and make sure the `FILTER` keyword (or
  `FITS_KEY_FILTER`) is set. `phcorr` uses it to check that flat and frame
  match.

### Synopsis

```bash
munipack flat [options] file(s)
```

### Options

#### Pre-corrections of the input flats

| Option | Argument | Default | Description |
|---|---|---|---|
| `-gain` | e⁻/ADU | header `GAIN`, otherwise 1 | Gain applied to every input flat. Overrides the header. |
| `-bias` | file | — | Master bias subtracted from every flat. |
| `-dark` | file | — | Master dark subtracted from every flat, scaled by *t*. |
| `-xdark` | factor | exposure-time ratio | Dark time factor *t*. If not given, *t* = `EXPTIME`(flat) / `EXPTIME`(dark), or 1 when unknown. |
| `-bitmask` | file | — | Mask frame: pixels ≤ 0.5 are not used and are set to 1 in the master. |
| `-st`, `--saturate` | ADU | header `SATURATE` | Upper limit of the raw values used (no gain applied). |
| `-th`, `--threshold` | ADU | machine epsilon | Lower limit of the raw values used (no gain applied). Must be > 0. |

#### Method and output

| Option | Argument | Default | Description |
|---|---|---|---|
| `--approximation` | `basic`, `standard` | `standard` | `standard`: the iterative photon-noise solution. `basic`: only the robust mean of the scaled frames, a rough estimate. |
| `-o`, `--output` | file | `flat.fits` | Output file. |
| `-B`, `--bitpix` | −32, 32, 16, 8 | −32 | Numerical type of the output. Use −32. Integer types scale the flat to a mean of 10², 10⁴ or 10⁹ (8, 16, 32 bits) and lose precision; `phcorr` normalises such flats. |
| `--verbose` | — | off | Prints the per-frame statistics and the gain estimate (see below). It also writes `GAIN_AVG` and `GAIN_STD` to the header. |
| `-h`, `--help` | — | — | Shows the help. |

### Gain estimate

For pure photon noise the corrected flat values satisfy *c*/σ² = 1 (in
electrons). With `--verbose`, `flat` prints for each frame the mean level
*c*, its standard error, the standard deviation σ and the ratio *g* =
*c*/σ², flagged as reliable or not. It also prints the mean of the ratios,
which is the **gain relative to the gain applied** (`-gain` or `GAIN`):

```text
Filename, mean level[ct], std.err., std.dev., gain, reliable:
flat-V_000006.fits:   99142.  0.151 313.73  1.007  T
...
 Estimated relative gain, std.err:  1.006  .002 (original gain was 2.300).
```

To estimate an unknown gain, start with `-gain 1`. Multiply the gain by the
reported ratio and repeat until the ratio is 1 within the accuracy you need.
Read noise, non-linearity and fixed-pattern structure all violate the
Poisson assumption, so check the result against another method, such as the
manufacturer's value or a photon-transfer curve.

### Outputs

- **Primary HDU:** the master flat, dimensionless, with a mean of 1.
- **`STDERR` extension:** the standard error of each pixel of the flat.
- **Header:**
  - `FILTER`: from the inputs.
  - `DATE-OBS`: of the first input.
  - `IMAGETYP`: copied from the inputs.
  - `GAIN_AVG`, `GAIN_STD`: the estimated relative gain and its spread, with
    `--verbose` only.
  - `COMMENT` lines: the inputs, the gain, bias, dark (with *t*) and mask
    used, `-th`, `-st` and the mean level.

### Examples

```bash
# V flat from bias- and dark-corrected twilight flats
munipack flat -bias bias.fits -dark dark30.fits -o flat_V.fits flat_V_*.fits

# gain given on the command line, saturation at 60000 ADU
munipack flat -gain 2.3 -dark d30.fits -st 6e4 -o f_R.fits f30_*R.fits

# estimate the gain (prints c, sigma and c/sigma^2 for every frame)
munipack flat --verbose -gain 1 -bias bias.fits -o flat_test.fits flat_V_*.fits

# quick, rough flat (robust mean only)
munipack flat --approximation basic -o flat_quick.fits flat_V_*.fits
```

(phcorr)=
## phcorr

Applies the photometric corrections (gain, bias, dark and flat-field) to
science frames and converts the values to electrons:

*C* = *g* (*S* − *t D* − *B*) / *F*

- *g* is the gain: from `-gain`, otherwise `GAIN`, otherwise 1 with a warning.
- *t* is the dark time factor: from `-xdark`, otherwise `EXPTIME`(frame) /
  `EXPTIME`(dark), otherwise 1.
- *F* is the flat-field, used as it is (a `flat` master has a mean of 1). It
  is divided by its robust mean with `--normalise-flat`, and always for
  integer flats.

A missing bias or dark counts as zero, and a missing flat as one. The masters
must have the dimensions of the frame. The flat must be positive everywhere.
`phcorr` warns when the flat has no filter. With `--verbose`, it also warns
when the filter of the flat differs from that of the frame, or when the
temperature of the dark differs by more than
`MUNIPACK_TEMPERATURE_TOLERANCE`.

The standard error of every corrected pixel goes to a `STDERR` extension:

σ² = [*g* max(*m*, 0) + *g*²(*r*² + *q* + τ²σ<sub>B</sub>² + *t*²σ<sub>D</sub>²)] / *F*² + (*C* σ<sub>F</sub> / *F*)²

- *m* is the detected signal *S* − *B*, from the mean of the eight
  neighbours of the pixel. It includes the sky, the stars and the dark
  current. The pixel itself is not used, so inverse-variance weights do not
  depend on its own noise.
- *r* is the read noise in ADU: `READNS`, or `RDNOISE` [e⁻] / *g*.
- *q* = 1/12 for integer frames (the quantisation noise).
- σ<sub>B</sub>, σ<sub>D</sub> and σ<sub>F</sub> are the `STDERR`s of the
  masters.
- τ = 1 − *t* when the dark was made with `dark -bias` and the same bias
  (`DARKBIAS`), otherwise 1.

The error is at least 0.01 *g*.

The output header carries the noise model in electrons, so `find`, `aphot`
and `psfphot` use the right one without further options.

### Synopsis

```bash
munipack phcorr [options] file(s)
```

### Options

#### Corrections

| Option | Argument | Default | Description |
|---|---|---|---|
| `-gain` | e⁻/ADU | header `GAIN`, otherwise 1 | Gain of the conversion to electrons. Overrides the header. Without either, the values stay in ADU (`BUNIT = 'ADU'`). |
| `-bias` | file | — | Master bias (from `bias`). |
| `-dark` | file | — | Master dark (from `dark`). |
| `-xdark` | factor | exposure-time ratio | Dark time factor *t*. Otherwise *t* = `EXPTIME`(frame) / `EXPTIME`(dark), or 1 when either is unknown. |
| `-flat` | file | — | Master flat-field (from `flat`). All values must be positive. |
| `--normalise-flat` | — | off | Divides the flat by its robust mean before use, so the frame keeps its mean level. Needed for a flat not normalised to 1, for example one from another program. Integer flats are always normalised. |

#### Masking

| Option | Argument | Default | Description |
|---|---|---|---|
| `-bitmask` | file | — | Mask frame of the same size: pixels of value 0 are replaced, all others (e.g. 1) are kept (`BITPIX = 8` is enough). Useful for hot pixels and bad columns. |
| `-xbitmask` | `ZERO`, `MEAN`, `MEDIAN` | `MEAN` | Replacement of masked pixels (case-insensitive; any other value acts as `ZERO`). `ZERO`: 0. `MEAN`: the robust mean of the corrected frame. `MEDIAN`: the robust mean of the unmasked pixels in a box around the pixel (falls back to `MEAN` when the box has none). |
| `-box` | pixels | 5 | Box size for `-xbitmask MEDIAN`. Use 3 for isolated masked pixels. |

#### Output files and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Writes the outputs to this directory (it must exist). The inputs are untouched. |
| `--enable-overwrite` | — | off | Allows replacing an existing output file, including the input itself (in-place correction). |
| `-O` | — | off | Turns on advanced output names (`--pattern`, `--mask`, `--format`); `--mask` alone also does. The names are written in the current directory: `-t` does not apply to them. |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against the input names; its groups are `\1`, `\2`, … (for `-O`). |
| `--mask` | string | — | Output name from the groups of `--pattern`, e.g. `'\1_c.\2'` (for `-O`). Unrelated to `-bitmask`. |
| `--format` | printf format | — | Numbered output names, counted from 0, e.g. `'cal_%03d.fits'` gives `cal_000.fits`, `cal_001.fits`, … (for `-O`). |
| `-B`, `--bitpix` | −32, 32, 16, 8 | −32 | Numerical type of the output. Integer types clip the values to 0 … 2<sup>BITPIX</sup> − 1. Keep −32 for frames in electrons. |
| `--verbose` | — | off | Prints the correction applied to each frame. |
| `-h`, `--help` | — | — | Shows the help. |

**Output file names.**

- **Default:** the input file name, without its directory, written in the
  current directory.
- **With `-t dir`:** `dir/` plus the input name.
- **With an existing output:** the run fails (exit status 1) unless
  `--enable-overwrite` is given. Running `phcorr` on files in the current
  directory therefore needs `-t`, `-O --mask …` or `--enable-overwrite`.
- **Lists:** `@list` with two names per line (`input output`) sets each
  output explicitly.

### Outputs

For every input frame, a FITS file with:

- **The first image HDU:** the corrected frame (in electrons when a gain is
  known). Its header is copied from the input and updated as follows:

  | Keyword | Value |
  |---|---|
  | `BUNIT` | `'COUNTS'` (electrons) when a gain is known, else `'ADU'` |
  | `GAIN` | 1: the values are electrons (when a gain is known) |
  | `PHCGAIN` | the gain used for the conversion [e⁻/ADU] (when a gain is known) |
  | `READNS`, `RDNOISE` | the read noise in electrons (if the input had either); without a gain only `READNS`, in ADU |
  | `SATURATE` | the saturation level converted to the output units |
  | `BIASLEV` | removed when a bias or dark is subtracted, otherwise converted |
  | `BITPIX` | the output type (`-B`) |
  | `COMMENT` | `PHCORR gain/flat-field/dark/bias/mask` records (the flat with its mean, the dark with *t*) |

- **`STDERR` extension:** the standard error of each pixel (see above). An
  existing `STDERR` of the input is replaced.
- **Other extensions** of the input (tables and so on) are copied unchanged.

### Examples

```bash
# full correction, outputs in cal/ with the same names
mkdir -p cal
munipack phcorr -bias bias.fits -dark dark300.fits -flat flat_V.fits \
    -t cal object_V_*.fits

# dark of the same exposure time (contains the bias), no bias frame
munipack phcorr -dark dark120.fits -flat flat_R.fits -t cal m51_R_*.fits

# outputs next to the inputs with a suffix: m51_001.fits -> m51_001_c.fits
munipack phcorr -bias bias.fits -dark dark.fits -flat flat_V.fits \
    -O --mask '\1_c.\2' m51_*.fits

# replace hot pixels by the local robust mean (3 x 3), correct in place
munipack phcorr -bias bias.fits -flat flat_V.fits -bitmask hot.fits \
    -xbitmask MEDIAN -box 3 --enable-overwrite frame_*.fits
```

```{note}
Older versions (before 0.6.2+hg1653 with these fixes) kept the raw `GAIN`,
`READNS` [ADU] and `BIASLEV` in the corrected header, and their `STDERR` had
a one-electron floor and no dark current. Frames corrected by them give
`find`, `aphot` and `psfphot` a wrong noise model. Correct the raw frames
again with this version.
```
