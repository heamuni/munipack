(phcorr)=
# phcorr

`phcorr` applies gain, bias, dark and flat corrections to science frames and
converts them to electrons.

Any master you leave out is simply skipped. The masters must match the frame
size, and the flat must be positive everywhere. The masters come from
{doc}`bias`, {doc}`dark` and {doc}`flat`.

## Usage

```bash
munipack phcorr [options] file(s)
```

## Options

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

## Output

One calibrated frame per input. Each gets a `STDERR` extension with the error
of every pixel, and its header carries the noise model in electrons. `find`,
`aphot` and `psfphot` pick this up automatically.

By default the output keeps the input's name and goes to the current
directory. `phcorr` refuses to overwrite an existing file, so when the inputs
are in the current directory use `-t`, `-O --mask …` or `--enable-overwrite`.

## Examples

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

(calibration-sequence)=
## Typical calibration sequence

The four routines turn raw CCD or CMOS exposures into calibrated frames in
electrons, ready for `find`, `aphot`, `psfphot` and `astrometry`.

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

## Notes

- Output names can also be set explicitly with an `@list` file of two names
  per line (`input output`).
- With `-O`, the outputs are written to the current directory; `-t` does not
  apply.
- Input lists, header keywords (`FITS_KEY_*`) and Docker use work as for the
  masters: see {ref}`masters-common`.
- Frames calibrated with older Munipack versions carry a wrong noise model.
  Calibrate the raw frames again with the current version before photometry.
