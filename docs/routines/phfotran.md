(phfotran)=

# phfotran

`phfotran` measures how your instrumental filters differ from a standard
photometric system and writes the transformation table used by
`phcal --tratab`.

Observe a standard field (e.g. Landolt or Stetson) with one frame per filter,
prepare the frames as for {doc}`phcal` (astrometry and photometry), and run
`phfotran` on them.

There are two models. The default `--model mag` fits the familiar colour
terms of the magnitudes, filter by filter, comparable to published
instrument colour terms. `--model flux` fits a linear mixing of photon
rates, which suits a filter that leaks light from a neighbouring band.

## Usage

```bash
munipack phfotran [options] file(s)
```

## Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | file | — | Reference catalogue of the standard field. |
| `-f`, `--filters` | list | — | Filters, one per frame (paired by `FILTER`). |
| `-C`, `--cal` | list | 1 | Factors multiplying the counts, in the order of `-f`. |
| `-E`, `--extin` | list | 0 | Extinction coefficients, for a transformation outside the atmosphere. |
| `--col-ra` | name | `RAJ2000` | Right Ascension column. |
| `--col-dec` | name | `DEJ2000` | Declination column. |
| `--col-mag` | list | — | Magnitude columns, in the order of `-f`. |
| `--col-magerr` | list | none | Magnitude error columns. |
| `--model` | `mag`, `flux` | `mag` | Colour terms of magnitudes, or linear mixing of photon rates. |
| `--color-order` | `1`, `2` | `1` | Order of the colour terms of the `mag` model. |
| `--photsys-ref` | name | — | Standard photometric system of the catalogue. |
| `--photsys-instr` | name | `PHOTSYS` of the frames | Name of your instrumental system. |
| `--area` | m² | `AREA` of the frames | Area of the telescope's entrance aperture. |
| `--phsystab` | file | bundled | Photometric system definitions. |
| `--list` | — | — | List the available photometric systems and exit. |
| `--tol` | deg | HWHM of the frames | Cross-match radius. |
| `--vocat` | file | bundled | VO catalogue configuration (column IDs). |
| `-o`, `--output` | file | `phfotran.fits` | Output table. |
| `--verbose` | — | off | Print the fitted coefficients and scatter. |
| `-h`, `--help` | — | — | Show the help. |

## Output

One FITS table (`phfotran.fits` unless `-o` says otherwise; an existing file
is overwritten) with the fitted transformation. Pass it to {doc}`phcal` with
`--tratab`.

## Examples

Fit the transformation on the standard field T Phe, then calibrate science
frames with it:

```bash
munipack phfotran --area 1.86 --photsys-instr DK154 --photsys-ref Johnson \
    -c T_Phe.fits --col-ra RA --col-dec DEC -f B,V,R,I --col-mag B,V,R,I \
    T_Phe_B.fits T_Phe_V.fits T_Phe_R.fits T_Phe_I.fits

munipack phcal --tratab phfotran.fits --photsys-instr DK154 --photsys-ref Johnson \
    -c field.fits --col-ra RA --col-dec DEC -f B,V,R,I --col-mag B,V,R,I \
    field_B.fits field_V.fits field_R.fits field_I.fits
```

Fit a photon-mixing matrix outside the atmosphere for a camera with a band
leak:

```bash
munipack phfotran --model flux -E 0.25,0.15,0.10 -c cat.fits --photsys-ref Johnson \
    -f B,V,R --col-mag Bmag,Vmag,Rmag -o tra_flux.fits B.fits V.fits R.fits
```

## Notes

- You need the same number of instrumental and standard filters.
- With `-E`, the airmass is computed from the observatory coordinates in the
  headers; without them it is undefined.
- For a strong band mixing, the `mag` model needs `--color-order 2`.
- The `--photsys-instr` name is stored in the table, and `phcal --tratab`
  checks that it matches your frames.
- Tables from older versions, written before the `mag` model existed, are
  read as `flux` tables.
