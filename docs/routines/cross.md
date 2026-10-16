(cross)=
# cross

`cross` cross-matches two FITS tables by position.

Each row of the second table gets the nearest row of the first, if it lies
within `--tol`. Use it to match two catalogues (e.g. from {doc}`cone`), or a
catalogue with the `PHOTOMETRY` table of a calibrated frame. The result is a
table of row indices, not a merged table. `cross` works offline.

## Usage

```bash
munipack cross --tol DEG [options] file1 file2
```

Each file is read from its first table extension. Choose another one with the
CFITSIO syntax, e.g. `'frame_cal.fits[PHOTOMETRY]'`.

## Options

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
of `cone` tables. Leave out `--col-pm-ra` and `--col-pm-dec`; {doc}`astrometry`
and {doc}`phcal` apply proper motions themselves.
```

## Output

The output table `CROSSMATCH` has two columns of 1-based row numbers: for
each row of file2, the matching row of file1, and vice versa (0 means no
match).

## Examples

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

## Notes

- The match is not strictly one-to-one: two close stars in one table can both
  pick the same star in the other.
