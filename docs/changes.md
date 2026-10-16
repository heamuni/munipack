# Differences from upstream

This fork is based on the upstream Munipack development version 0.6.2+hg1653 (Mercurial revision 1653). It fixes a long list of bugs and adds a few new tools and options. This page gives the short version. The full technical list, with how each bug was found and checked, is in [BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md), and how the methods were tested is in {doc}`development`.

Most of the photometry fixes matter on faint skies of a few electrons per pixel. On typical ground-based frames, results change by less than their noise.

## What was fixed

- **`artificial`:** a realistic detector noise model (Poisson sky, read noise, gain, bias); frames usable again with `--scale` in pixels per degree; correct star positions, PSF centring, sky gradient and airmass; reproducible output with `--seed`.
- **`find`:** the threshold now means what it says on Poisson noise; fewer false detections near bright stars; faint stars no longer rejected by the shape tests.
- **`aphot`, `find`, `psfphot`:** the sky is no longer underestimated on faint frames, so faint stars no longer come out too bright.
- **`aphot`:** correct errors (the gain is used), a true `HWHM` keyword, and apertures contaminated by a bright neighbour are flagged.
- **`gphot`:** rebuilt on the `aphot` enclosed energy, with correct errors and sky; `--legacy` keeps the old method.
- **`phcal`:** correct physical quantities (FLAM, FNU, ABMAG, STMAG ...); frames matched to filters by `FILTER`; no spurious colour correction; `-C` and `-q PHOTRATE` work.
- **`phfotran`:** now actually runs.
- **`phcorr`:** `-box` is honoured, and errors and output headers are correct.
- **`dark`:** faint master darks no longer come out low.
- **`kombine`:** correct resampling grid and edges, `--pm-ra`/`--pm-dec` work, `STDERR` on the right scale.
- **`astrometry`:** mirrored and rotated frames solve; distortion fits correctly; wrong matches on dense or saturated frames are rejected; frames with too few stars fail on their own without stopping the run.
- **`fits`:** `--update` or `--remove-keys` on an extension no longer deletes the file.
- **Command line:** exit status is 0 on success and 1 on error, with the reason on stderr; unknown or incomplete options give an error instead of being silently ignored.

## What is new

- **`psfphot`:** PSF photometry with star groups and an empirical PSF.
- **`aphot --apcorr`:** aperture correction to total flux, also when the PSF changes over the field.
- **`phcal`:** a systematic error term (`--syserr auto`) and colour terms (`--color-term`).
- **`phfotran --model mag`:** a magnitude transformation model.
- **`find --false-rate`:** set the detection threshold as a false-alarm rate.
- **`astrometry`:** matching from the header WCS or a known pointing, also with an unknown roll; a scale prior for the blind match; parallel processing.
- **`kombine --flux-conserve`:** flux-conserving resampling.
- **`artificial`:** `--dark`, `--readnoise`, `--bias`, `--seed`, `--no-watermark`.

## Compatibility notes

These changes give different results than older versions. The complete list is in [BUGS.md §5](https://github.com/heamuni/munipack/blob/develop/BUGS.md#5-compatibility-notes-behaviour-that-changed).

```{warning}
Re-run `aphot` on frames processed by an older version before using `phcal`. Old tables have a wrong `HWHM` and wrong errors, which make `phcal` pick the wrong aperture.
```

- **`--scale`** is in pixels per degree in `artificial`, `astrometry` and `kombine` (an upstream change).
- **`artificial`** frames now have real noise, so they differ from older versions. Note that `--hwhm` is σ for the GAUSS and SEEING PSFs and α for MOFFAT, not a half width (as upstream).
- **`find`** finds fewer noise peaks on faint skies; re-tune `-th` if you rely on a given depth.
- **`phcal`** adds a systematic error by default (`--syserr 0` gives the old errors) and uses the aperture correction when present (`--apcorr none` ignores it).
- **`phfotran`** fits magnitudes by default; `--model flux` gives the old tables.
- **`astrometry`** needs at least five stars, rejects insignificant matches, and tries the header WCS or pointing before a blind match (`--blind` forces the blind match). It also uses the header scale as a prior; `--scale-free` turns that off.
- **`kombine` and `dark`** use a clipped mean by default; `--hampel` gives the old robust mean.
- **`phcorr`** output in electrons has `GAIN = 1`, with the original gain kept in `PHCGAIN`.
- **`psfphot`** flags stars with S/N below 2 (FLAG 256).
- **Exit status** is 0 or 1, so scripts that tested other codes need updating.
