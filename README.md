# README
## Article Information
This repository provides access to the data and source code used for the manuscript  
**Intentional or mistaken? Termites start pairing broadly but refine choices**  
**Nobuaki Mizumoto, Elijah P. Carroll**  
**Paper DOI:** [TBA]

This study investigates the broadness of tandem running behavior in a termite _Reticulitermes amamianus_. During mating season, mating pairs perform tandem running while looking for a nest site, by males following females. We investigated whether males follow similar to females but non-mating individuals (workers and soldiers). Behavioral observations include posture tracking of laboratory recordings of behavioral interactions. This repository includes data and the Python/R scripts.

## Table of Contents
This repository includes tracking data, R code to analyze it, and Python code for video analysis. Trajectory data are not included due to file size.

* [README](./README.md)
* [code](./code)
  * [`data_processing_well.py`](./code/data_processing_well.py), [`data_processing_dish.py`](./code/data_processing_dish.py) - Format `.h5` files (well experiments) to `.feather` with interpolation and smoothing, for well and dish experiment respectively.
  * [`helper_function.py`](./code/helper_function.py) - functions used for above data_processig
  * [`processing.R`](./code/processing.R) - Format trajectories for visualization and statistical analysis
  * [`output.R`](./code/output.R) - Visualization and statistics for experiments
  * [`source.R`](./code/source.R) - for loading packages and functions

## Setup & Dependencies
Scripts are written in R and Python, tested on Windows 11 (64-bit).

### R Session Info
```
```

### Python Environment
```
```

## Citation
TBA

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Nobuaki Mizumoto: nzm0095@auburn.edu
