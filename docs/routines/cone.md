(cone)=
# cone

`cone` downloads the stars around a position from a VizieR catalogue.

It asks a VizieR server for all stars within a radius of a position and saves
them as a FITS table or a VOTable, sorted by magnitude. This is how you get
the reference catalogue that {doc}`astrometry`, {doc}`phcal` and
{doc}`artificial` need. Pass the catalogue's column names to those routines
(for UCAC4: `RAJ2000`, `DEJ2000`, `pmRA`, `pmDE`, `Vmag`).

`cone` needs network access. To get the coordinates of an object by its name,
use {doc}`sesame`; to convert a VOTable or draw a chart, use {doc}`votable`;
to match two catalogues, use {doc}`cross`.

## Usage

```bash
munipack cone [options] [--] RA DEC
munipack cone --list-catalogues | --list-servers
```

RA and Dec are in decimal degrees. Put `--` before them when the declination
is negative, so that it is not read as an option.

## Options

**Catalogue and server**

| Option | Default | Description |
|---|---|---|
| `-c`, `--cat` *alias* | `UCAC4` | Predefined catalogue; a prefix is enough (`-c Gaia`). |
| `--id` *ID* | — | Any VizieR catalogue by its identifier, e.g. `V/154`. |
| `--server` *alias* | CDS | VizieR mirror (see `--list-servers`). |
| `--url` *URL* | — | Full query URL instead of catalogue and server. |
| `--vocat` *file* | installed `VOcat_conf.xml` | Catalogue definitions (aliases, IDs, column names, mirrors). |
| `--list-catalogues` | — | List the predefined catalogues and their columns. |
| `--list-servers` | — | List the VizieR mirrors. |

**Selection**

| Option | Default | Description |
|---|---|---|
| `-r`, `--radius` *deg* | 0.1 (6 arcmin) | Search radius. |
| `--magmin` *mag* | — | Bright limit: keep only fainter stars. |
| `--magmax` *mag* | — | Faint limit: keep only brighter stars. |
| `-s`, `--sort` *column* | catalogue magnitude | Sort by this column; the magnitude limits apply to it. |
| `--par` *key=value* | — | Extra VizieR constraint, e.g. `'-out.max=500'`. |
| `--Johnson-patch` | off | UCAC4 only: convert Gunn r, i to Johnson R, I. |

**Output**

| Option | Default | Description |
|---|---|---|
| `-o`, `--output` *file* | `cone.fits` | Output file; a suffix other than `.fits` gives a VOTable. |
| `--type` *fits\|votable* | by suffix | Output format. |
| `--verbose` | off | Print the query and the server's answer. |
| `-h`, `--help` | — | Show the help. |

## Output

The FITS output is a single binary table with every catalogue column plus
`_r`, the distance from the centre. The epoch of the positions is stored in
the `JEPOCH` header keyword when the catalogue provides it. A VOTable output
can be converted later with {doc}`votable`.

## Examples

UCAC4 stars within 0.3° of M67:

```bash
munipack cone -r 0.3 -o m67.fits -- 132.85 11.81
```

Gaia DR3 stars brighter than G = 17 at a southern position:

```bash
munipack cone -c Gaia --magmax 17 -r 0.2 -o field.fits -- 245.9 -26.5
```

Any VizieR catalogue by its ID. Munipack does not know its magnitude column,
so name it with `--sort` for the limit to work:

```bash
munipack cone --id V/154 -s umag --magmax 20 -r 0.2 -o sdss.fits -- 180.0 0.0
```

## Notes

- The aliases of `-c` and `--server` are case-sensitive. If both `-c` and
  `--id` are given, `-c` wins.
- Only the last `--par` is used; join several constraints with `&`.
- To give your own catalogue an alias, copy `VOcat_conf.xml`, add a row to its
  `CATALOGUE` table and pass the copy with `--vocat`.
- When nothing matches, `cone` says so and writes no file.
