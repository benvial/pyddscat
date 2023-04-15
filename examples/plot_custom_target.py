#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/pyddscat


"""
Custom targets
=====================

Define your targets from functions.
"""

# sphinx_gallery_thumbnail_number = 1

import numpy as np
import pyddscat as pd
import matplotlib.pyplot as plt


##############################################################################
# Let's build a torus


R = 1
r = R / 5


def torus(x, y, z, R=R, r=r):
    """Build a torus"""
    rsq = (R - np.sqrt(x**2 + y**2)) ** 2 + z**2
    return np.where(rsq < r**2, 1, 0)


D = R + 2 * r
target = pd.targets.FROM_FILE.fromfunction(
    torus, (-D, -D, -D), (D, D, D), d=r / 4, material=2 - 0.1j
)


##############################################################################
# Plot target
target.show()
