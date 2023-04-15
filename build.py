#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/pyddscat

import setuptools
from distutils import log

import sys
import subprocess


class FortranBuild(setuptools.Command):
    """Customized setuptools build command."""

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        log.info("Building fortran...")
        command = ["make", "all"]
        if subprocess.call(command) != 0:
            sys.exit(-1)
        super().run()
