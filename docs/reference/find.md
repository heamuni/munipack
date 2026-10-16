# find

`find` detects stars on calibrated frames. Run it after calibration and
before photometry or astrometry: it stores the star list in a `FIND` table
inside each frame, which `aphot`, `psfphot` and `astrometry` then read.

The method follows DAOFIND from Stetson's DAOPHOT. The frame is smoothed with
a Gaussian of the stars' FWHM, and peaks above the threshold become
candidates. Shape tests (sharpness and roundness) then drop hot pixels,
cosmic rays and bad columns, and a candidate in the halo of a brighter star
must stand out above that halo.

## Synopsis

```bash
munipack find [options] file(s)
```

Files can also be given as a list, `@list.txt`, or as `-` to read names from
standard input.

## Options

### Detection

| Option | Default | Description |
|---|---|---|
| `-f`, `--fwhm` | 3 | FWHM of the stars in pixels; a much smaller value splits stars and finds noise. |
| `-th`, `--threshold` | 7 | Detection threshold in sigmas (see {ref}`find-threshold`). |
| `--false-rate` | — | Set the threshold by the number of noise detections you accept per frame; overrides `-th`. |
| `--lothresh` | 5 | Ignore pixels more than this many sigmas below the sky. |
| `--saturate` | header, else from the data type | Saturation level in ADU; saturated areas are not reported as stars. |
| `--read-noise` | from header | Read noise in ADU. |

### Shape tests

| Option | Default | Description |
|---|---|---|
| `--shape-sigma` | 3 | Accept sharpness and roundness within this many sigmas of a typical star; `0` uses the fixed limits only. |
| `--shrplo` | 0.2 | Lower sharpness limit (fixed). |
| `--shrphi` | 1.0 | Upper sharpness limit (fixed). |
| `--rndlo` | −1.0 | Lower roundness limit (fixed). |
| `--rndhi` | 1.0 | Upper roundness limit (fixed). |

### Output and general

| Option | Description |
|---|---|
| `-t`, `--target-directory` | Write the results to this directory and leave the input frames untouched. |
| `-O` | Turn on advanced output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` | Regular expression matched against input names, default `(.+)\.(.+)` (for `-O`). |
| `--mask` | Output name built from the pattern groups, e.g. `'\1_f.\2'` (for `-O`). |
| `--format` | Numbered output names, e.g. `'out_%03d.fits'` (for `-O`). |
| `--verbose` | Print the sky level, the number of stars and the threshold for each frame. |
| `-h`, `--help` | Show the help. |

(find-threshold)=
## Choosing the threshold

The threshold sets how many noise peaks sneak in as "stars". Every pixel is
tested, so even perfect noise gives some false detections. The default
`-th 7` is safe and picks only fairly bright stars; `-th 5` goes deeper and
leaves about one noise peak on a 2048 × 2048 frame; `-th 4` reaches fainter
stars but lets in on the order of a hundred.

If you would rather say what you can tolerate, use `--false-rate`. For
example, `--false-rate 1` picks the threshold that gives about one noise
detection per frame, whatever the frame size.

When the header has `GAIN`, `find` accounts for photon noise on a faint sky,
so the threshold means the same on a dark space-based frame as on a bright
ground-based one.

```{warning}
A stacked frame that still carries the `GAIN` of a single exposure gets a
threshold that is too high. Correct `GAIN` before running `find`.
```

The shape tests do not remove noise peaks, because near the threshold they
look just like faint stars, so the threshold alone sets the contamination.

## Shape tests

By default a candidate passes if its sharpness and roundness lie within the
fixed limits, or within `--shape-sigma` sigmas of a typical star on the same
frame. The second rule keeps faint stars that photon noise has pushed outside
DAOPHOT's classic limits. To use the fixed limits only, give `--shape-sigma 0`
(both tests), `--shrplo`/`--shrphi` (sharpness) or `--rndlo`/`--rndhi`
(roundness).

## Output

`find` adds (or replaces) a `FIND` table extension in each frame; the image
and its header are left as they are. Without `-t` or `-O` the input files are
updated in place.

The table lists the stars brightest first, with columns `X`, `Y` (position in
pixels), `PEAKRATIO` (peak height relative to the threshold), `SHARP` and
`ROUND`. Its header records the settings used, such as `FWHM`, `THRESH` and
`SATURATE`.

## Examples

Detect stars with the default threshold:

```bash
munipack find -f 3.5 frame.fits
```

A deeper search allowing about one noise detection per frame, written to
`found/` so the originals stay untouched:

```bash
mkdir -p found
munipack find -f 4.9 --false-rate 1 --verbose -t found frame_*.fits
```

Give the saturation and read noise when the header lacks them, then measure
the stars:

```bash
munipack find -f 3 -th 5 --saturate 60000 --read-noise 3.7 frame.fits
munipack aphot frame.fits
```

## Notes

Frames made by `artificial` have a small watermark in the top-right corner,
which `find` picks up as a few "stars". Simulate with `--no-watermark` to
avoid it.
