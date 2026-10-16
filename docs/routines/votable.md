(votable)=
# votable

`votable` converts a VOTable to a FITS table, CSV, text or an SVG finding chart.

The input is a VOTable from `cone -o file.xml` (see {doc}`cone`) or from any
VO tool; it must be a VOTable, not FITS. `votable` works offline.

## Usage

```bash
munipack votable [options] -o output file.xml
```

## Options

| Option | Default | Description |
|---|---|---|
| `-o`, `--output` *file* | `votable.fits` | Output file; always give it. |
| `--type` *fits\|csv\|txt\|svg\|xml* | by suffix | Output format. |
| `--verbose` | off | Verbose messages. |
| `-h`, `--help` | — | Show the help. |

**SVG chart**

| Option | Default | Description |
|---|---|---|
| `-pt` *projection* | `TAN` | Projection (only gnomonic is implemented). |
| `-pa` *deg* | — | RA of the chart centre. |
| `-pd` *deg* | — | Dec of the chart centre. |
| `-ps` *px/deg* | — | Chart scale. |
| `-pw` *px* | 500 | Canvas width (use together with `-ph`). |
| `-ph` *px* | 500 | Canvas height (use together with `-pw`). |
| `-ml` *mag* | — | Reference magnitude for symbol size; fainter stars are drawn smaller. |
| `--col-ra` *field* | — | RA field, e.g. `RAJ2000`. |
| `--col-dec` *field* | — | Dec field, e.g. `DEJ2000`. |
| `--col-mag` *field* | — | Magnitude field, e.g. `f.mag`. |

## Output

One file in the format chosen by `--type` or by the suffix of `-o`: a FITS
table, CSV, tab-separated text, an SVG chart, or a VOTable again.

## Examples

Convert a cone search to a FITS table or a CSV file:

```bash
munipack votable -o pleiades.fits pleiades.xml
munipack votable -o pleiades.csv pleiades.xml
```

An 800 × 600 px finding chart of M67:

```bash
munipack cone -r 0.1 -o m67.xml -- 132.85 11.81
munipack votable -o m67.svg -pa 132.85 -pd 11.81 -ps 4000 -ml 16 -pw 800 -ph 600 \
    --col-ra RAJ2000 --col-dec DEJ2000 --col-mag f.mag m67.xml
```

## Notes

- For a chart, always give `-pa`, `-pd`, `-ps`, `-ml` and the three `--col-*`
  names: they have no usable defaults, and the field names must match the
  VOTable exactly.
- North is up and east is left on the chart.
