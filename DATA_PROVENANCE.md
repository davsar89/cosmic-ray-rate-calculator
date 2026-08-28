# Data provenance, attribution, and limitations

## Source graphic

The original project states that `data.csv` was digitized from the Wikimedia Commons graphic **“Cosmic ray flux versus particle energy.”** The authoritative file page is:

```text
https://commons.wikimedia.org/wiki/File:Cosmic_ray_flux_versus_particle_energy.svg
```

Attribution recorded on that page:

- title: *Cosmic ray flux versus particle energy*;
- author: Sven Lafebre;
- date: 8 January 2019 for the current 2018-version graphic;
- source: the author’s own work, after S. Swordy (2001) and A. De Angelis and M. Pimenta (2018).

The source page offers several licenses. For the local rasterized copy and any copyrightable adaptation of it in this repository, use the **Creative Commons Attribution-ShareAlike 3.0 Unported** option unless the maintainers deliberately select another offered license:

```text
https://creativecommons.org/licenses/by-sa/3.0/
```

That license requires attribution, a link to the license, an indication of changes, and share-alike licensing for adaptations. The Wikimedia file page remains authoritative for the complete authorship history and all available license choices.

## Local changes and derived files

The included file `Cosmic_ray_flux_versus_particle_energy.svg.png` is a `1920 x 2304` PNG rasterization of the source SVG. The original archive did not document the rasterization tool or any changes beyond the format/rendering conversion.

`data.csv` is a numerical digitization of the plotted spectrum. This corrective patch:

- preserved all 74 supplied numeric data rows;
- replaced the comment-only pseudo-header with a defined CSV header;
- added comments explaining provenance and precision;
- changed no energy or flux value.

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

## Checksums of the original supplied files

These hashes identify the unmodified files in the uploaded archive before this corrective patch changed the CSV header and comments:

```text
SHA-256  data.csv
3ddf276dcdadb912310607d4e11f9a60a5a2cbc068f22ec6db1692459b77271b

SHA-256  Cosmic_ray_flux_versus_particle_energy.svg.png
82cffef47be9bc4afa2aa9c18fa982db03f034604bac39456fd1edf917be1aeb
```

The patched loader also computes a semantic SHA-256 identifier from the sorted binary floating-point values. That identifier is unaffected by comments or row order and prevents an approximation from being silently used with different numeric data.

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

A scientifically stronger replacement should cite a primary source, preserve reported uncertainties, describe any unit conversions, and record a reproducible transformation from the source data to this CSV schema.
