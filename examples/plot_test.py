#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat


"""
A simple example
================

Homogeneous, isotropic finite cylinder with hemispherical endcaps.
"""


import os
import numpy as np
import pyddscat as pd
import matplotlib.pyplot as plt


##############################################################################
# Parameters

length = 0.100
radius = 0.020

material = os.path.join(os.path.abspath(os.path.curdir), "Au_evap")
target = pd.targets.CYLNDRCAP(length, radius, d=0.004, material=material)

##############################################################################
# Plot target
target.show()

##############################################################################
# Create a job to be run
job = pd.DDscat(target=target)

##############################################################################
# Change the range of calculated wavelengths and ambient index
job.settings.wavelengths = pd.ranges.How_Range(0.300, 0.600, 31)
job.settings.NAMBIENT = 1.0

##############################################################################
# Run the job
job.calculate()

out = job.output

##############################################################################
# Check optical theorem
optical_theorem = out["Q_abs"] + out["Q_sca"] - out["Q_ext"]
print(optical_theorem)
assert np.allclose(optical_theorem, 0, atol=1e-4)


##############################################################################
# Plot
ax = out.plot(["Q_sca", "Q_abs", "Q_ext"])
ax.legend(loc=0)
ax.set_xlabel("wavelength (microns)")
ax.set_ylabel("Scattering cross sections")
