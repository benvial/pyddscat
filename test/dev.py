#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat


def test_show():
    import pyddscat as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import os.path

    # Establish target geometry (in um)
    length = 0.100
    radius = 0.020

    def torus(x, y, z, R=0.5, r=0.1):
        """Build a torus"""
        rsq = (R - np.sqrt(x**2 + y**2)) ** 2 + z**2
        return np.where(rsq < r**2, 1, 0)

    target = pd.targets.FROM_FILE.fromfunction(torus, (-1, -1, -1), (1, 1, 1))

    target.VTRconvert()
    # target.show(screenshot="test.png",off_screen=True)
    target.show(off_screen=True)
