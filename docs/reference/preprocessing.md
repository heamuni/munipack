# bias, dark, flat, phcorr — calibration frames and corrections

These four actions turn raw CCD or CMOS exposures into calibrated frames in
electrons, ready for `find`, `aphot`, `psfphot` and `astrometry`.

| Action | What it makes |
|---|---|
| `bias` | a master bias from zero-exposure frames |
| `dark` | a master dark from dark frames |
| `flat` | a master flat (mean of 1) for one filter |
| `phcorr` | calibrated science frames, one output per input |

(calibration-sequence)=
## Typical calibration sequence

```bash
# 1. master bias (more than 10 frames is better)
munipack bias -o bias.fits bias_*.fits

# 2. master dark with the bias removed, so it can be scaled to any exposure time
munipack dark -bias bias.fits -o dark120.fits dark120_*.fits

# 3. one master flat per filter
munipack flat -bias bias.fits -dark dark120.fits -o flat_V.fits flat_V_*.fits

# 4. calibrate the science frames into cal/, leaving the originals untouched
mkdir -p cal
munipack phcorr -bias bias.fits -dark dark120.fits -flat flat_V.fits \
    -t cal object_V_*.fits
```

For every pixel, `phcorr` computes *C* = *g* (*S* − *t D* − *B*) / *F*: the
raw frame *S* minus the bias *B* and the dark *D* (scaled by the exposure-time
ratio *t*), divided by the flat *F* and multiplied by the gain *g*.

You can work in one of two ways:

- **Darks of the same exposure time as your science frames, no bias.** Run
  `dark` and `phcorr` without `-bias`. The dark then includes the bias. This
  is the simplest choice when all frames share one exposure time.
- **Bias plus scaled darks.** Make the dark with `-bias` and give the same
  bias to `flat` and `phcorr`. One long dark then serves every exposure time.

```{warning}
Do not scale a dark that still contains the bias: the bias would be scaled
too. Use such a dark only for frames of its own exposure time.
```

## Common to bias, dark and flat

**Pixel rejection.** Raw pixel values outside the range set by `-th`
(threshold) and `-st` (saturation) are left out of the average. Both are in
raw ADU. Saturation defaults to the `SATURATE` header keyword.

**Masks.** `-bitmask` takes a FITS image of the same size. Pixels with value
1 are used and pixels with value 0 are skipped.

**Header checks.** `bias` and `dark` warn when the inputs differ in exposure
time or chip temperature, and `flat` warns when the filters differ. If your
files use other keyword names, set the `FITS_KEY_*` environment variables,
for example `FITS_KEY_EXPTIME=EXPOSURE`. Others are `FITS_KEY_FILTER`,
`FITS_KEY_TEMPERATURE`, `FITS_KEY_DATEOBS`, `FITS_KEY_GAIN` and
`FITS_KEY_SATURATE`. The read noise is taken from `READNS` (in ADU) or
`RDNOISE` (in electrons).

**Input lists.** Give files as arguments, as `@list.txt` (one name per line),
or as `-` to read the names from standard input.

**Output.** Each master has a `STDERR` extension with the error of each
pixel. Keep the default float output (`-B -32`); an existing output file is
replaced.

(bias)=
## bias

Averages zero-exposure frames into a master bias, the electronic offset of the
readout. Each pixel is a robust mean, so cosmic-ray hits are rejected. Bias
frames are cheap, so take many.

```bash
munipack bias [options] file(s)
```

| Option | Default | Description |
|---|---|---|
| `-o`, `--output` | `bias.fits` | Output file. |
| `-st`, `--saturate` | header `SATURATE` | Upper limit of raw values used [ADU]. |
| `-th`, `--threshold` | 1 | Lower limit of raw values used [ADU]. |
| `-bitmask` | — | Mask frame (0 = skip pixel). |
| `-B`, `--bitpix` | −32 | Output type: −32, 32, 16 or 8. |
| `--verbose` | off | Print the mean level, exposure time and temperature. |
| `-h`, `--help` | — | Show the help. |

```bash
# master bias from all bias frames in the directory
munipack bias -o bias.fits bias_*.fits

# frames listed in a file, with a bad-column mask
munipack bias -bitmask badcols.fits -o bias.fits @bias.lst
```

(dark)=
## dark

Averages dark frames into a master dark. With `-bias`, the bias is subtracted
first, so the master holds the dark current alone and `flat` and `phcorr` can
scale it to other exposure times. Remember that dark current depends strongly
on chip temperature, so take darks at the temperature of your science frames.

By default each pixel is a mean clipped by the expected noise, which rejects
cosmic rays and hot outliers while staying accurate for very weak darks.

```bash
munipack dark [-bias file] [--arith | --hampel] [options] file(s)
```

| Option | Default | Description |
|---|---|---|
| `-bias` | — | Master bias subtracted from every dark frame. |
| `--arith` | off | Plain arithmetic mean, no rejection. |
| `--hampel` | off | Robust (Hampel) mean, the old default; underestimates weak darks. Not with `--arith`. |
| `-o`, `--output` | `dark.fits` | Output file. |
| `-st`, `--saturate` | header `SATURATE` | Upper limit of raw values used [ADU]. |
| `-th`, `--threshold` | 1 | Lower limit of raw values used [ADU]. |
| `-bitmask` | — | Mask frame (0 = skip pixel). |
| `-B`, `--bitpix` | −32 | Output type; keep −32 for a bias-subtracted dark, whose values can be below zero. |
| `--verbose` | off | Print the noise model, mean level, exposure time and temperature. |
| `-h`, `--help` | — | Show the help. |

```bash
# dark of the science exposure time, bias left in
munipack dark -o dark120.fits dark120_*.fits

# bias-subtracted dark, to be scaled to any exposure time
munipack dark -bias bias.fits -o dark300.fits dark300_*.fits
```

(flat)=
## flat

Builds a master flat from twilight, dome or night-sky flats. The flat maps the
response of the whole light path: pixel sensitivity, vignetting, dust and
filter non-uniformity. Each input is corrected for bias and dark, then the
frames are combined using their photon statistics, which suppresses cosmic
rays and other defects. The result is normalised to a mean of 1.

For a good flat, expose to about half of the full well, take more than 13
frames, and make one flat per filter with the `FILTER` keyword set, so that
`phcorr` can check that flat and frame match.

```bash
munipack flat [options] file(s)
```

| Option | Default | Description |
|---|---|---|
| `-gain` | header `GAIN`, else 1 | Gain [e⁻/ADU]; overrides the header. |
| `-bias` | — | Master bias subtracted from every flat. |
| `-dark` | — | Master dark subtracted from every flat. |
| `-xdark` | exposure-time ratio | Dark scaling factor *t*. |
| `-bitmask` | — | Mask frame (0 = skip pixel; set to 1 in the master). |
| `-st`, `--saturate` | header `SATURATE` | Upper limit of raw values used [ADU]. |
| `-th`, `--threshold` | tiny positive value | Lower limit of raw values used [ADU]. |
| `--approximation` | `standard` | `standard` (full fit) or `basic` (quick robust mean). |
| `-o`, `--output` | `flat.fits` | Output file. |
| `-B`, `--bitpix` | −32 | Output type; keep −32, integer types lose precision. |
| `--verbose` | off | Print per-frame statistics and a gain estimate. |
| `-h`, `--help` | — | Show the help. |

```bash
# V flat from bias- and dark-corrected twilight flats
munipack flat -bias bias.fits -dark dark30.fits -o flat_V.fits flat_V_*.fits

# estimate the gain of the camera
munipack flat --verbose -gain 1 -bias bias.fits -o flat_test.fits flat_V_*.fits
```

**Estimating the gain.** With `--verbose`, `flat` reports the gain relative to
the one it applied. Start with `-gain 1`, multiply by the reported ratio and
repeat until the ratio is close to 1. Read noise and non-linearity bias this
estimate, so compare it with the manufacturer's value.

(phcorr)=
## phcorr

Applies gain, bias, dark and flat corrections to science frames and converts
them to electrons. Any master you leave out is simply skipped. The masters
must match the frame size, and the flat must be positive everywhere.

Each output gets a `STDERR` extension with the error of every pixel, and its
header carries the noise model in electrons. `find`, `aphot` and `psfphot`
pick this up automatically.

```bash
munipack phcorr [options] file(s)
```

### Corrections

| Option | Default | Description |
|---|---|---|
| `-gain` | header `GAIN`, else 1 | Gain [e⁻/ADU]; without one, values stay in ADU. |
| `-bias` | — | Master bias. |
| `-dark` | — | Master dark. |
| `-xdark` | exposure-time ratio | Dark scaling factor *t* (1 when exposure times are unknown). |
| `-flat` | — | Master flat. |
| `--normalise-flat` | off | Divide the flat by its mean first; use for flats not normalised to 1. |

### Masking

| Option | Default | Description |
|---|---|---|
| `-bitmask` | — | Mask frame; pixels with value 0 are replaced (e.g. hot pixels, bad columns). |
| `-xbitmask` | `MEAN` | Replacement value: `ZERO`, `MEAN` (frame mean) or `MEDIAN` (local mean in a box). |
| `-box` | 5 | Box size for `MEDIAN`; use 3 for isolated pixels. |

### Output

| Option | Default | Description |
|---|---|---|
| `-t`, `--target-directory` | — | Write outputs to this (existing) directory. |
| `--enable-overwrite` | off | Allow replacing existing files, including the inputs. |
| `-O` | off | Enable custom output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` | `(.+)\.(.+)` | Regular expression matched against input names (for `-O`). |
| `--mask` | — | Output name built from the pattern groups, e.g. `'\1_c.\2'` (for `-O`). |
| `--format` | — | Numbered output names, e.g. `'cal_%03d.fits'` (for `-O`). |
| `-B`, `--bitpix` | −32 | Output type; keep −32 for frames in electrons. |
| `--verbose` | off | Print the corrections applied to each frame. |
| `-h`, `--help` | — | Show the help. |

By default the output keeps the input's name and goes to the current
directory. `phcorr` refuses to overwrite an existing file, so when the inputs
are in the current directory use `-t`, `-O --mask …` or `--enable-overwrite`.

```bash
# full correction, outputs in cal/ with the same names
munipack phcorr -bias bias.fits -dark dark300.fits -flat flat_V.fits \
    -t cal object_V_*.fits

# outputs next to the inputs with a suffix: m51_001.fits -> m51_001_c.fits
munipack phcorr -bias bias.fits -dark dark.fits -flat flat_V.fits \
    -O --mask '\1_c.\2' m51_*.fits

# replace hot pixels by the local mean, correcting in place
munipack phcorr -bias bias.fits -flat flat_V.fits -bitmask hot.fits \
    -xbitmask MEDIAN -box 3 --enable-overwrite frame_*.fits
```

## Notes

- `flat` keeps all input frames in memory (about 8 bytes per pixel per
  frame), so a large set of big frames needs plenty of RAM.
- Output names can also be set explicitly with an `@list` file of two names
  per line (`input output`).
- With `-O`, the outputs are written to the current directory; `-t` does not
  apply.
- In Docker, mount your data on `/data` and pass environment variables with
  `-e`, e.g. `-e FITS_KEY_EXPTIME=EXPOSURE`. Add `-i` to `docker run` when
  reading file names from standard input.
- Frames calibrated with older Munipack versions carry a wrong noise model.
  Calibrate the raw frames again with the current version before photometry.
