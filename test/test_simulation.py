#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat

import pyddscat as pd
import matplotlib.pyplot as plt
import numpy as np
import os.path


def test_simulation():
    # Establish target geometry (in um)
    length = 0.100
    radius = 0.020

    material = os.path.join(os.path.dirname(__file__), "Au_evap")
    target = pd.targets.CYLNDRCAP(length, radius, d=0.004, material=material)

    # Create a job to be run
    job = pd.DDscat(target=target)

    # Change the range of calculated wavelengths and ambient index
    job.settings.wavelengths = pd.ranges.How_Range(0.300, 0.600, 2)
    job.settings.NAMBIENT = 1.0

    # Run the job locally
    job.calculate(silent=True)

    out = job.output
    optical_theorem = out["Q_abs"] + out["Q_sca"] - out["Q_ext"]
    print(optical_theorem)
    assert np.allclose(optical_theorem, 0, atol=1e-4)

    # ax = out.plot(["Q_sca", "Q_abs", "Q_ext"])
    # ax.legend(loc=0)
    # ax.set_xlabel("wavelength (microns)")
    # ax.set_ylabel("Scattering cross sections")
    # # plt.show()
