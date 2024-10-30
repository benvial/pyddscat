#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat


"""
Scattering by a cylinder
========================

Homogeneous, isotropic finite cylinder with hemispherical endcaps.
"""


import numpy as np
import pyddscat as pd
import matplotlib.pyplot as plt


##############################################################################
# Parameters
wl = 0.5
a = 5 * wl / (2 * np.pi)
d = a / 2.81
material = 0.96 + 1.01j
target = pd.targets.Sphere(a, d=d, material=material)

##############################################################################
# Plot target
target.show()

##############################################################################
# Create a job to be run
job = pd.DDscat(target=target)

##############################################################################
# Change the range of calculated wavelengths and ambient index
job.settings.wavelengths = pd.ranges.How_Range(wl, wl, 1)
job.settings.NRFLD = True  # Make nearfield calc
job.settings.scat_planes = [pd.ranges.Scat_Range(0, 0, 180, 1)]
# job.settings.NRFLD_EXT = np.array(
#     (0.5,) * 6
# )  # Set nearfield box to 0.5x on all sides

##############################################################################
# Run the job
job.calculate()

##############################################################################
# Compute fields
lines = [pd.fileio.Line((-1.5 * a, 0, 0), (1.5 * a, 0, 0), 501)]
job.postprocess_fields(lines=lines)


x = job.fields["x"]

normE2 = (
    np.abs(job.fields["Ex"]) ** 2
    + np.abs(job.fields["Ey"]) ** 2
    + np.abs(job.fields["Ez"]) ** 2
)

plt.plot(x / a, normE2)
plt.show()
# job.fields
