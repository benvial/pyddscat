#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat


"""
Custom targets
=====================

Define your targets from functions.
"""

# sphinx_gallery_thumbnail_number = 1

import numpy as np
import pyddscat as pd

##############################################################################
# Let's build a torus


def torus(x, y, z, R=0.5, r=0.1):
    """Build a torus"""
    rsq = (R - np.sqrt(x**2 + y**2)) ** 2 + z**2
    return np.where(rsq < r**2, 1, 0)


target = pd.targets.FROM_FILE.fromfunction(torus, (-1, -1, -1), (1, 1, 1))


##############################################################################
# Plot with pyvista

target.show()
