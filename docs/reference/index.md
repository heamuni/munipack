# Reference

Every Munipack task is an action, run as `munipack <action> [options] file(s)`.
Each page below covers one action or a few related ones: what it does, its
options and a couple of examples. For the full option list on your own
machine, run `munipack <action> --help`.

| Page | Actions | What for |
|---|---|---|
| {doc}`preprocessing` | `bias`, `dark`, `flat`, `phcorr` | Make calibration frames and apply them |
| {doc}`find` | `find` | Detect stars |
| {doc}`aphot` | `aphot` | Aperture photometry |
| {doc}`psfphot` | `psfphot` | PSF photometry |
| {doc}`gphot` | `gphot` | Growth-curve photometry |
| {doc}`astrometry` | `astrometry` | Astrometric calibration (WCS) |
| {doc}`phcal` | `phcal`, `phfotran` | Photometric calibration and colour transformations |
| {doc}`kombine` | `kombine`, `colouring` | Stacks, mosaics and colour images |
| {doc}`artificial` | `artificial` | Simulated frames |
| {doc}`catalogues` | `cone`, `sesame`, `votable`, `cross` | Catalogues and the Virtual Observatory |
| {doc}`utilities` | `fits`, `timeseries` | FITS files and time series |

`munipack` on its own takes three options: `--version` prints the version and
licence, `-h`/`--help` lists the actions, and `--verbose` makes the action report
its progress. An action exits with status 0 on success and 1 on an error,
with the reason printed to stderr, so it is easy to use in scripts. If you
run Munipack from the container, see {doc}`../install`.

```{toctree}
:hidden:

preprocessing
find
aphot
psfphot
gphot
astrometry
phcal
kombine
artificial
catalogues
utilities
```
