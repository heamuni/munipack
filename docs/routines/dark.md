(dark)=
# dark

`dark` averages dark frames into a master dark.

With `-bias`, the bias is subtracted first, so the master holds the dark
current alone and `flat` and `phcorr` can scale it to other exposure times
(see {ref}`calibration-sequence`). Dark current depends strongly on chip
temperature, so take darks at the temperature of your science frames.

By default each pixel is a mean clipped by the expected noise, which rejects
cosmic rays and hot outliers while staying accurate for very weak darks.

## Usage

```bash
munipack dark [-bias file] [--arith | --hampel] [options] file(s)
```

## Options

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

## Output

The master dark, with a `STDERR` extension holding the error of each pixel.
Give it to `flat` and `phcorr` with `-dark`.

## Examples

```bash
# dark of the science exposure time, bias left in
munipack dark -o dark120.fits dark120_*.fits

# bias-subtracted dark, to be scaled to any exposure time
munipack dark -bias bias.fits -o dark300.fits dark300_*.fits
```

## Notes

- Do not scale a dark that still contains the bias: the bias would be scaled
  too. Use such a dark only for frames of its own exposure time.
- Pixel rejection, masks, header checks and input lists work as for all
  masters: see {ref}`masters-common`.
