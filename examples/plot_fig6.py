#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat


"""
Scattering by a sphere
======================

Reproduces the top panel of Fig. 6 in DDSCAT manual.
"""


import numpy as np
import pyddscat as pd
import matplotlib.pyplot as plt


##############################################################################
# Parameters

a = 1
d = a / 10
material = 1.7 + 0.1j
target = pd.targets.Sphere(a, d=d, material=material)

##############################################################################
# Plot target
# target.show()

##############################################################################
# Create a job to be run
job = pd.DDscat(target=target)

##############################################################################
# Change the range of calculated wavelengths and ambient index
job.settings.wavelengths = pd.ranges.How_Range(
    2 * np.pi * a / 8, 2 * np.pi * a / 0.01, 100, "INV"
)
job.settings.scat_planes = [pd.ranges.Scat_Range(0, 0, 0, 1)]

##############################################################################
# Run the job
job.calculate()

##############################################################################
# Plot
out = job.output
wls = job.settings.wavelengths.table
x = 2 * np.pi * a / wls

fig, ax = plt.subplots(figsize=(10, 4))
plt.plot(x, out["Q_sca"], label="Q_sca")
plt.plot(x, out["Q_abs"], label="Q_abs")
plt.plot(x, out["Q_ext"], label="Q_ext")
ax.legend(loc=0)
ax.set_xlabel(r"$x = 2\pi a/\lambda$")
plt.yscale("log")
plt.xlim(0, 8)
plt.ylim(1e-2)
plt.tight_layout()

plt.show()
