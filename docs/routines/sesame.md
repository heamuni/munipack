(sesame)=
# sesame

`sesame` turns an object name into RA and Dec.

It looks up the name with the CDS Sesame service (SIMBAD, NED, VizieR) and
prints the coordinates in degrees. Use it to feed {doc}`cone` or
`astrometry --rcen/--dcen` (see {doc}`astrometry`) without typing coordinates.
`sesame` needs network access.

## Usage

```bash
munipack sesame [options] "object name"
```

## Options

| Option | Default | Description |
|---|---|---|
| `object` | — | Object name (required); quote names with spaces. |
| `--verbose` | off | Print the request and the full answer. |
| `-h`, `--help` | — | Show the help. |

## Output

Two numbers on standard output: RA and Dec in decimal degrees. No file is
written.

## Examples

```bash
munipack sesame "beta per"
# 47.0422186 40.9556467
```

Resolve a name and search around it in one go:

```bash
munipack cone -r 0.25 -o ngc5139.fits -- $(munipack sesame "NGC 5139")
```
