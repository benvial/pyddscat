

<!-- <a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/-/releases" target="_blank"><img alt="Release" src="https://img.shields.io/endpoint?url=https://gitlab.com/phokaia/pyddscat/-/jobs/artifacts/main/raw/logobadge.json?job=badge&labelColor=c9c9c9"></a>  -->
<!-- <a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/commits/main" target="_blank"><img alt="Release" src="https://img.shields.io/gitlab/pipeline/phokaia/pyddscat/main?logo=gitlab&labelColor=dedede&style=for-the-badge"></a>  -->
<!-- <a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/commits/main" target="_blank"><img alt="Release" src="https://img.shields.io/gitlab/coverage/phokaia/pyddscat/main?logo=python&logoColor=e9d672&style=flat-square"></a> -->
<a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/commits/main" target="_blank"><img alt="pipeline" src="https://gitlab.com/phokaia/pyddscat/badges/main/pipeline.svg?style=flat-square"></a> 
<a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/commits/main" target="_blank"><img alt="coverage" src="https://gitlab.com/phokaia/pyddscat/badges/main/coverage.svg?style=flat-square"></a> 
<a class="reference external image-reference" href="https://black.readthedocs.io/en/stable/" target="_blank"><img alt="Release" src="https://img.shields.io/badge/code%20style-black-dedede.svg?logo=python&logoColor=e9d672&style=flat-square"></a>
<a class="reference external image-reference" href="https://gitlab.com/phokaia/pyddscat/-/blob/main/LICENSE.txt" target="_blank"><img alt="Release" src="https://img.shields.io/badge/license-GPLv3-blue?color=f8742d&logo=open-access&logoColor=f8742d&style=flat-square"></a>



# pyddscat -- A Python interface to DDSCAT

`pyddscat` is a Python package for interfacing to the popular scattering simulator
[DDSCAT](http://www.astro.princeton.edu/~draine/DDSCAT.html). `pyddscat` provides a rich toolset to:

* Create standard `DDSCAT` scattering targets based on physical (rather than dipole dimensions)
* Construct and visualize complex custom scattering targets
* Manage the job parameters found in the `ddscat.par` file
* Organize iterative jobs requiring multiple targets or input parameters
* Script job submission to cluster queue managers
* Maintain profiles and defaults for deployment on platforms other than the local machine
* Load, plot and manipulate `DDSCAT` output tables
* Manage the output from multiple jobs through results collections
* Work with and visualize near field results as multidimensional `numpy` arrays
* Suitable for interactive or scripted use



## Installation

You will need a `fortran` compiler, e.g.:

```
mamba install -c conda-forge fortran-compiler
```

Then clone and install locally with `pip`:

```
git clone git@gitlab.com:phokaia/pyddscat.git
cd pyddscat
pip install -e .
```

## Documentation

See [phokaia.gitlab.io/pyddscat](https://phokaia.gitlab.io/pyddscat).

## Examples

See the [example gallery](https://phokaia.gitlab.io/pyddscat/examples).

## Acknowledgement

This project is based on [ScatPy](https://pythonhosted.org/ScatPy/index.html) and has been ported to Python 3.

