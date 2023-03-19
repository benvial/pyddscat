#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat

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
    job.settings.wavelengths = pd.ranges.How_Range(0.300, 0.600, 15)
    job.settings.NAMBIENT = 1.0

    # Run the job locally
    job.calculate(silent=True)

    out = job.output

    optical_theorem = out["Q_abs"] + out["Q_sca"] - out["Q_ext"]
    print(optical_theorem)
    assert np.allclose(optical_theorem, 0, atol=1e-4)

    ax = out.plot(["Q_sca", "Q_abs", "Q_ext"])
    ax.legend(loc=0)
    ax.set_xlabel("wavelength (microns)")
    ax.set_ylabel("Scattering cross sections")
    # plt.show()


def test_show():
    def torus(x, y, z, R=0.5, r=0.1):
        """Build a torus"""
        rsq = (R - np.sqrt(x**2 + y**2)) ** 2 + z**2
        return np.where(rsq < r**2, 1, 0)

    target = pd.targets.FROM_FILE.fromfunction(torus, (-1, -1, -1), (1, 1, 1))
    # target.show(screenshot="test.png",off_screen=True)
    target.show(off_screen=True)
