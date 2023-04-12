

# pyddscat -- A Python interface to DDSCAT

`pyddscat` is a Python package for interfacing to the popular scattering simulator
[DDSCAT](http://www.astro.princeton.edu/~draine/DDSCAT.html). pyddscat provides a rich toolset to:

* Create standard `DDSCAT` scattering targets based on physical (rather than dipole) dimensions
* Construct and visualize complex custom scattering targets
* Manage the job parameters found in the `ddscat.par` file
* Organize iterative jobs requiring multiple targets or input parameters
* Script job submission to cluster queue managers
* Maintain profiles and defaults for deployment on platforms other than the local machine
* Load, plot and manipulate `DDSCAT` output tables
* Manage the output from multiple jobs through results collections
* Work with and visualize nearfield results as multidimensional `numpy` arrays
* Suitable for interactive or scripted use



## Installation

You will need a fortran compiler, e.g.:

```
mamba install -c conda-forge fortran-compiler
```

Then install locally with `pip`:

```
cd path/to/pyddscatt
pip install -e .
```

## Documentation

## Download


## Examples

## Acknowledgenent

This project is based on [ScatPy](https://pythonhosted.org/ScatPy/index.html) and has been ported to Python 3.

