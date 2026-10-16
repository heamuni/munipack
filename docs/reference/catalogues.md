# Catalogues and Virtual Observatory

Four small actions fetch reference catalogues and work with them:

| Action | What it does |
|---|---|
| [`cone`](#cone) | downloads the stars around a position from a VizieR catalogue |
| [`sesame`](#sesame) | turns an object name into RA and Dec |
| [`votable`](#votable) | converts a VOTable to FITS, CSV, text or an SVG chart |
| [`cross`](#cross) | cross-matches two FITS tables by position |

`cone` and `sesame` need network access; `votable` and `cross` work offline.

(cone)=
## cone

`cone` asks a VizieR server for all stars within a radius of a position and
saves them as a FITS table or a VOTable, sorted by magnitude. This is how you get the reference catalogue
that `astrometry`, `phcal` and `artificial` need. Pass the catalogue's column
names to those actions (for UCAC4: `RAJ2000`, `DEJ2000`, `pmRA`, `pmDE`,
`Vmag`).

```bash
munipack cone [options] [--] RA DEC
munipack cone --list-catalogues | --list-servers
```

RA and Dec are in decimal degrees. Put `--` before them when the declination
is negative, so that it is not read as an option.

### Options

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

The FITS output is a single binary table with every catalogue column plus
`_r`, the distance from the centre. The epoch of the positions is stored in
the `JEPOCH` header keyword when the catalogue provides it.

### Examples

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

### Notes

- The aliases of `-c` and `--server` are case-sensitive. If both `-c` and `--id` are
  given, `-c` wins.
- Only the last `--par` is used; join several constraints with `&`.
- To give your own catalogue an alias, copy `VOcat_conf.xml`, add a row to its
  `CATALOGUE` table and pass the copy with `--vocat`.
- When nothing matches, `cone` says so and writes no file.

(sesame)=
## sesame

`sesame` looks up an object name with the CDS Sesame service (SIMBAD, NED,
VizieR) and prints its RA and Dec in degrees. Use it to feed `cone` or
`astrometry --rcen/--dcen` without typing coordinates.

```bash
munipack sesame [options] "object name"
```

| Option | Default | Description |
|---|---|---|
| `object` | — | Object name (required); quote names with spaces. |
| `--verbose` | off | Print the request and the full answer. |
| `-h`, `--help` | — | Show the help. |

### Examples

```bash
munipack sesame "beta per"
# 47.0422186 40.9556467
```

Resolve a name and search around it in one go:

```bash
munipack cone -r 0.25 -o ngc5139.fits -- $(munipack sesame "NGC 5139")
```

(votable)=
## votable

`votable` converts a VOTable (from `cone -o file.xml` or any VO tool) to a FITS
table, CSV, tab-separated text or an SVG finding chart. The input must be a
VOTable, not FITS.

```bash
munipack votable [options] -o output file.xml
```

### Options

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

For a chart, always give `-pa`, `-pd`, `-ps`, `-ml` and the three `--col-*`
names: they have no usable defaults, and the field names must match the
VOTable exactly. North is up and east is left.

### Examples

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

(cross)=
## cross

`cross` pairs the rows of two FITS tables by position: each row of the second
table gets the nearest row of the first, if it lies within `--tol`. Use it to
match two catalogues, or a catalogue with the `PHOTOMETRY` table of a
calibrated frame. The result is a table of row indices, not a merged table.

```bash
munipack cross --tol DEG [options] file1 file2
```

Each file is read from its first table extension. Choose another one with the
CFITSIO syntax, e.g. `'frame_cal.fits[PHOTOMETRY]'`.

### Options

| Option | Default | Description |
|---|---|---|
| `--tol` *deg* | near zero | Match radius; always give it (the default matches nothing), e.g. `3e-4` (about 1″). |
| `--col-ra` *col1,col2* | `RA,RA` | RA columns of file1 and file2; give both names. |
| `--col-dec` *col1,col2* | `DEC,DEC` | Dec columns of file1 and file2; give both names. |
| `--col-pm-ra` *col1,col2* | — | Proper motion in RA (mas/yr); do not use, see below. |
| `--col-pm-dec` *col1,col2* | — | Proper motion in Dec (mas/yr); do not use, see below. |
| `--col-epoch-ra` *column* | — | Accepted but ignored. |
| `--col-epoch-dec` *column* | — | Accepted but ignored. |
| `--col-mag` *col1,col2* | — | Magnitude columns; not used by the match. |
| `--ftol` *ratio* | 1 | Flux tolerance; not used by the match. |
| `--vocat` *file* | — | Accepted but ignored. |
| `-o`, `--output` *file* | `cross.fits` | Output file (overwritten). |
| `--verbose` | off | Verbose messages. |
| `-h`, `--help` | — | Show the help. |

```{warning}
The proper-motion correction is broken in this version and spoils the match
of `cone` tables. Leave out `--col-pm-ra` and `--col-pm-dec`; `astrometry` and
`phcal` apply proper motions themselves.
```

The output table `CROSSMATCH` has two columns of 1-based row numbers: for
each row of file2, the matching row of file1, and vice versa (0 means no
match). The match is not strictly one-to-one: two close stars in one table
can both pick the same star in the other.

### Examples

UCAC4 against Gaia DR3 over the same field:

```bash
munipack cross --tol 3e-4 --col-ra RAJ2000,RA_ICRS --col-dec DEJ2000,DE_ICRS \
    -o ucac4_gaia.fits ucac4.fits gaia.fits
```

The stars of a calibrated frame against the reference catalogue:

```bash
munipack cross --tol 3e-4 --col-ra RAJ2000,RAJ2000 --col-dec DEJ2000,DEJ2000 \
    -o matched.fits 'frame_cal.fits[PHOTOMETRY]' ucac4.fits
```
