(bias)=
# bias

`bias` averages zero-exposure frames into a master bias, the electronic offset
of the readout.

It is the first step of the {ref}`calibration sequence <calibration-sequence>`.
Each pixel is a robust mean, so cosmic-ray hits are rejected. Bias frames are
cheap, so take many (more than 10 is better).

## Usage

```bash
munipack bias [options] file(s)
```

## Options

| Option | Default | Description |
|---|---|---|
| `-o`, `--output` | `bias.fits` | Output file. |
| `-st`, `--saturate` | header `SATURATE` | Upper limit of raw values used [ADU]. |
| `-th`, `--threshold` | 1 | Lower limit of raw values used [ADU]. |
| `-bitmask` | — | Mask frame (0 = skip pixel). |
| `-B`, `--bitpix` | −32 | Output type: −32, 32, 16 or 8. |
| `--verbose` | off | Print the mean level, exposure time and temperature. |
| `-h`, `--help` | — | Show the help. |

## Output

The master bias, with a `STDERR` extension holding the error of each pixel.
Give it to `dark`, `flat` and `phcorr` with `-bias`.

## Examples

```bash
# master bias from all bias frames in the directory
munipack bias -o bias.fits bias_*.fits

# frames listed in a file, with a bad-column mask
munipack bias -bitmask badcols.fits -o bias.fits @bias.lst
```

(masters-common)=
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

**Docker.** Mount your data on `/data` and pass environment variables with
`-e`, e.g. `-e FITS_KEY_EXPTIME=EXPOSURE`. Add `-i` to `docker run` when
reading file names from standard input.
