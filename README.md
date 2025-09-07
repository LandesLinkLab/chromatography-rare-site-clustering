# Rare adsorption site clustering in chromatography — MATLAB code

This repository contains MATLAB code used to simulate how **rare adsorption site clustering** affects peak broadening in chromatography.

> Primary reference:  Nikita Kovalenko and Christy F. Landes, The Impact of Rare Adsorption Site Clustering on Peak Broadening in Chromatography, Analytical Chemistry (2025). https://doi.org/10.1021/acs.analchem.5c03108

## Repository structure

```
src/                    % main functions
  single_run.m          % core simulation for one set of inputs
  script_column_run.m   % driver script

```

## Quick start

```matlab
run src/script_column_run.m
```
This will produce a histogram of retention times for a small synthetic setup.


## License

MIT (see `LICENSE`).

## Citation

Please cite the paper:
Nikita Kovalenko and Christy F. Landes, The Impact of Rare Adsorption Site Clustering on Peak Broadening in Chromatography, Analytical Chemistry (2025). https://doi.org/10.1021/acs.analchem.5c03108
