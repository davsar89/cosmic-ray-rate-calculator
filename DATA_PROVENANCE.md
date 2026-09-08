# Data provenance, attribution, and limitations

## Source graphic

The original project states that `data.csv` was digitized from Sven Lafebre's Wikimedia Commons graphic [Cosmic ray flux versus particle energy](https://commons.wikimedia.org/wiki/File:Cosmic_ray_flux_versus_particle_energy.svg).

Attribution recorded on that page:

- title: *Cosmic ray flux versus particle energy*;
- author: Sven Lafebre;
- date: 8 January 2019 for the current 2018-version graphic;
- source: the author’s own work, after S. Swordy (2001) and A. De Angelis and M. Pimenta (2018).

The source page offers several licenses. For the local rasterized copy and any copyrightable adaptation of it in this repository, the selected option is [Creative Commons Attribution-ShareAlike 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/).

That license requires attribution, a link to the license, an indication of changes, and share-alike licensing for adaptations. The Wikimedia file page remains authoritative for the complete authorship history and all available license choices.

## Local changes and derived files

The included file `Cosmic_ray_flux_versus_particle_energy.svg.png` is a `1920 x 2304` PNG rasterization of the source SVG. The original archive did not document the rasterization tool or any changes beyond the format/rendering conversion.

`data.csv` preserves all 74 originally supplied numeric data rows. Only the header and provenance comments have been clarified; no energy or flux value has changed.

To the extent that `data.csv`, `rate_comparison.png`, or other generated visual material is legally an adaptation of the source graphic, redistribute that material under CC BY-SA 3.0 with the attribution above. This statement is a provenance and reuse notice, not legal advice.

## What is not available

The supplied archive did not contain:

- a primary experimental data table;
- a digitization script or project file;
- original pixel/source coordinates for the extracted points;
- a record of digitizer settings or transformations;
- statistical or systematic uncertainties;
- a stable identifier for the exact source revision used beyond the included raster file.

Because of these omissions, the numeric values should be treated as approximate readings from a plot. Their stored decimal length reflects digitizer output, not physical measurement accuracy.

The source labels flux per GeV and particle energy in eV but supplies no species decomposition or explicit kinetic-versus-total energy convention. The calculator uses those plotted coordinates without composition or per-nucleon transformations. See [Scientific scope and limitations](README.md#scientific-scope-and-limitations) for the physical interpretation and independent physics references.

An approximation retains its calibration spectrum. Comparisons check the numeric spectrum and tail convention directly, so row order and comments do not affect matching. File edits are read on the next call, including for the default approximation.

## Replacing the spectrum

A replacement CSV must use this header:

```csv
energy_eV,flux_per_m2_sr_s_GeV
```

Each row must contain a finite, strictly positive energy and flux. Row order is arbitrary; duplicate energies are rejected.

When a spectrum is modified, regenerate the project artifacts with:

```bash
python generate_artifacts.py
```

For command-line comparisons against a custom dataset, explicitly request recalibration:

```bash
python cosmic_ray_rate.py 1e15 \
  --csv replacement.csv \
  --recalibrate-approximation \
  --check-formula
```

If the replacement spectrum does not span the built-in approximation break points, recalibration restricts them to the dataset's energy range and prints a note; explicit break energies can be supplied with `--approx-breaks-ev`.

A scientifically stronger replacement should cite a primary source, preserve reported uncertainties, describe any unit conversions, and record a reproducible transformation from the source data to this CSV schema.
