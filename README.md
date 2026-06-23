# README
## Article Information
This repository provides access to the source code used for the manuscript  
**Beyond mistakes: same-sex partner acceptance and mating filter coexist in termite pairing**  
**Nobuaki Mizumoto, Elijah P. Carroll**  
**Paper DOI:** [TBA]

This study investigates the broadness of tandem running behavior in a termite _Reticulitermes amamianus_. During mating season, mating pairs perform tandem running while looking for a nest site, by males following females. We investigated whether males follow similar to females but non-mating individuals (workers and soldiers). Behavioral observations include posture tracking of laboratory recordings of behavioral interactions. This repository includes the Python/R scripts.

## Table of Contents
This repository includes R and Python code for data analysis. Data will be available elsewhere.

* [README](./README.md)
* [code](./code)
  * [`data_processing_well.py`](./code/data_processing_well.py), [`data_processing_dish.py`](./code/data_processing_dish.py) - Format `.h5` files (well experiments) to `.feather` with interpolation and smoothing, for well and dish experiments, respectively.
  * [`helper_function.py`](./code/helper_function.py) - functions used for above data_processig
  * [`processing.R`](./code/processing.R) - Format trajectories for visualization and statistical analysis
  * [`output.R`](./code/output.R) - Visualization and statistics for experiments
  * [`source.R`](./code/source.R) - for loading packages and functions

## Setup & Dependencies
Scripts are written in R and Python, tested on Windows 11 (64-bit).

### R Environment

The R environment is managed using `renv`. Package versions used in the analyses are recorded in `renv.lock`.
To reproduce the R environment:

```r
install.packages("renv")
renv::restore()
```

### Python Environment

Python scripts were developed using Python 3.11.4
Required Python packages can be installed using:

```python
pip install -r requirements.txt
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Nobuaki Mizumoto: nzm0095@auburn.edu
