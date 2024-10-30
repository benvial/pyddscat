#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat

import pyddscat as pd
import numpy as np


def test_show():
    def torus(x, y, z, R=0.5, r=0.1):
        """Build a torus"""
        rsq = (R - np.sqrt(x**2 + y**2)) ** 2 + z**2
        return np.where(rsq < r**2, 1, 0)

    target = pd.targets.FROM_FILE.fromfunction(torus, (-1, -1, -1), (1, 1, 1))
    # target.show(screenshot="test.png",off_screen=True)
    target.show(off_screen=True)
