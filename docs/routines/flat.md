(flat)=
# flat

`flat` builds a master flat, normalised to a mean of 1, from twilight, dome or
night-sky flats of one filter.

The flat maps the response of the whole light path: pixel sensitivity,
vignetting, dust and filter non-uniformity. Each input is corrected for bias
and dark, then the frames are combined using their photon statistics, which
suppresses cosmic rays and other defects.

For a good flat, expose to about half of the full well, take more than 13
frames, and make one flat per filter with the `FILTER` keyword set, so that
`phcorr` can check that flat and frame match.

## Usage

```bash
munipack flat [options] file(s)
```

## Options

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

## Output

The master flat (mean 1), with a `STDERR` extension holding the error of each
pixel. Give it to `phcorr` with `-flat`.

## Examples

```bash
# V flat from bias- and dark-corrected twilight flats
munipack flat -bias bias.fits -dark dark30.fits -o flat_V.fits flat_V_*.fits

# estimate the gain of the camera
munipack flat --verbose -gain 1 -bias bias.fits -o flat_test.fits flat_V_*.fits
```

## Estimating the gain

With `--verbose`, `flat` reports the gain relative to the one it applied.
Start with `-gain 1`, multiply by the reported ratio and repeat until the
ratio is close to 1. Read noise and non-linearity bias this estimate, so
compare it with the manufacturer's value.

## Notes

- `flat` keeps all input frames in memory (about 8 bytes per pixel per
  frame), so a large set of big frames needs plenty of RAM.
- Pixel rejection, masks, header checks and input lists work as for all
  masters: see {ref}`masters-common`.
