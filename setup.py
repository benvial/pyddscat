#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat

import setuptools

import setuptools
from distutils import log

import sys
import subprocess


from setuptools.command.build_ext import build_ext


class FortranBuild(build_ext):
    """Customized setuptools build command."""

    def run(self):
        log.info("Building fortran...")
        command = ["make", "all"]
        if subprocess.call(command) != 0:
            sys.exit(-1)
        super().run()


setuptools.setup(
    cmdclass={
        "build_ext": FortranBuild,
    },
    scripts=[
        "bin/ddscat",
        "bin/vtrconvert",
        "bin/ddpostprocess",
        "bin/calltarget",
    ],
)
