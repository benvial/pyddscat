
Documentation
-------------



    Version |release| |  built on |today|




Introduction
============

pyddscat is a Python package for interfacing to the popular scattering simulator
`DDSCAT <https://ddscat.wikidot.com/>`_. pyddscat provides a rich toolset to:

* Create standard DDSCAT scattering targets based on physical (rather than dipole) dimensions
* Construct and visualize complex custom scattering targets
* Manage the job parameters found in the ddscat.par file
* Organize iterative jobs requiring multiple targets or input parameters
* Script job submission to cluster queue managers
* Maintain profiles and defaults for deployment on platforms other than the local machine
* Load, plot and manipulate DDSCAT output tables
* Manage the output from multiple jobs through results collections
* Work with and visualize nearfield results as multidimensional numpy arrays
* Suitable for interactive or scripted use

Installation
============


1. Start by getting a working version of `Python <http://www.python.org/getit/>`_. 
   We strongly suggest to use `Anaconda <https://anaconda.org/>`_.

2. Clone the gitlab repository:

.. code:: 

    git clone git@gitlab.com:phokaia/pyddscat.git

and install the package locally

.. code:: 

    cd pyddscat
    pip install -e .

This should compile DDSCAT binaries and install the requirements for pyddscat. You should have a working fortran compiler.

Getting Started
===============
Once you have ``pyddscat`` installed, start by reading the :doc:`userguide/userguide`
and reviewing the :doc:`examples`. More advanced users will benefit from referring
to the :doc:`api/api`. 

Acknowledgements and Citing
===========================

If this code contributes to a publication, please cite

    
    The Phokaia Developers; pyddscat: Python package for interacting with DDSCAT.
    Available at https://gitlab.com/phokaia/pyddscat.


License
=======
pyddscat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see `<http://www.gnu.org/licenses/>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

