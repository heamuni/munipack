# colouring

`colouring` makes a colour image from frames taken through standard filters,
typically Johnson B, V and R.

It transforms the frames into a colour space such as CIE 1931 XYZ, so the
result approximates natural colours. The output is a 3-D FITS cube with one
plane per band.

## Usage

```bash
munipack colouring [options] files
```

List the frames from short to long wavelength (B, V, R).

All bands must be on the same pixel grid. Resample them with {doc}`kombine`
first, using the same `--rcen`, `--dcen`, `--width`, `--height`, `--scale`,
`--angle` and `--reflex` for each band.

For correct colours the bands need a white balance. Calibrated frames, or
ctph values given with `-q`, give absolute intensities. Otherwise give
weights (`-w`), a white spot or a white star; without any of them
`colouring` warns that the colours are not balanced.

## Options

Lists are comma-separated, with one item per input file (use semicolons if
your locale uses a decimal comma).

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cspace-input` | name | guessed from FILTER | Input colour space, e.g. `'Johnson BVR'`. |
| `--cspace-output` | name | `'CIE 1931 XYZ'` | Output colour space. |
| `--list` | — | — | List the available colour-space transformations and exit. |
| `--ctable` | file | installed `ctable.dat` | Colour-transformation definitions. |
| `--phsystab` | file | installed `photosystems.fits` | Photometric systems table. |
| `-w`, `--weights` | `w1,w2,…` | — | Weights of the bands. |
| `-q`, `--ctphs` | `ctph1,ctph2,…` | from calibrated frames | ctph of each band, for calibrated output. |
| `-b`, `--backs` | `b1,b2,…` | estimated | Background level of each band. |
| `--disable-back` | — | off | Do not estimate or subtract backgrounds. |
| `--white-spot` | `x,y` [pix] | — | Centre of a region that should come out white. |
| `--white-star` | `x,y` [pix] | — | Centre of a star taken as white (B − V = 0). |
| `--white-radius` | real [pix] | 7 | Radius of the white spot or star. |
| `-o`, `--output` | file | `colour.fits` | Output file (overwritten if it exists). |
| `-B`, `--bitpix` | `-32`, `32`, `16`, `8` | `-32` | Data type of the output. |
| `--verbose` | — | off | Print the guessed colour space and progress. |
| `-h`, `--help` | — | — | Show the help and exit. |

## Output

One 3-D FITS cube with one plane per band.

## Examples

Turn B, V, R frames on a common grid into a colour image:

```bash
munipack colouring -o m27.fits m27_B.fits m27_V.fits m27_R.fits
```

Balance uncalibrated frames on a white star at pixel (512, 488):

```bash
munipack colouring --white-star 512,488 --white-radius 5 -o m27.fits \
    m27_B.fits m27_V.fits m27_R.fits
```

See which colour-space transformations are available:

```bash
munipack colouring --list
```

## Notes

- The input colour space is guessed only for three frames with FILTER
  exactly `B`, `V`, `R` (or `X`, `Y`, `Z`); other sets need `-c`.
