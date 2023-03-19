#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

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
)
# if __name__ == "__main__":
#     setuptools.setup()
#     #     cmdclass={
#     #         "build_ext": Build,
#     #     },
#     #     # scripts=["bin/ddscat", "bin/ddpostprocess", "bin/vtrconvert", "bin/calltarget"],
#     # )
