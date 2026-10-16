# Reference

One page per action or group of related actions. Each lists the options with
their units and defaults (as printed by `munipack <action> --help`), the
method, the outputs and examples.

All actions are run as `munipack <action> [options] file(s)`; in the
container, `munipack` is the entry point (see {doc}`../install`).

| Page | Actions | Purpose |
|---|---|---|
| {doc}`preprocessing` | `bias`, `dark`, `flat`, `phcorr` | Master calibration frames and their application |
| {doc}`find` | `find` | Star detection |
| {doc}`aphot` | `aphot` | Aperture photometry and aperture corrections |
| {doc}`psfphot` | `psfphot` | PSF photometry |
| {doc}`gphot` | `gphot` | Growth-curve photometry |
| {doc}`astrometry` | `astrometry` | Astrometric calibration (WCS) |
| {doc}`phcal` | `phcal`, `phfotran` | Photometric calibration and instrument transformations |
| {doc}`kombine` | `kombine`, `colouring` | Stacks, mosaics and colour frames |
| {doc}`artificial` | `artificial` | Simulated frames |
| {doc}`catalogues` | `cone`, `sesame`, `votable`, `cross` | Catalogues and Virtual Observatory |
| {doc}`utilities` | `fits`, `timeseries` | FITS files and time series |

The front end itself:

```text
Usage: munipack [--version] [-h] [--verbose] [action]
```

`--version` prints the version and licence, `--help` the list of actions,
`--verbose` makes the actions log their progress. The exit status is 0 on
success and 1 on any error, with the reason on stderr.

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
