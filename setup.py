#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/pyddscat


import shutil
import sys
import subprocess
from contextlib import suppress
from pathlib import Path
from setuptools import Command, setup
from setuptools.command.build import build
from setuptools import logging


class CustomCommand(Command):
    def initialize_options(self) -> None:
        self.bdist_dir = None
        self.pkg_name = None

    def finalize_options(self) -> None:
        self.pkg_name = self.distribution.get_name().replace("-", "_")
        with suppress(Exception):
            self.bdist_dir = Path(self.get_finalized_command("bdist_wheel").bdist_dir)

    def run(self) -> None:
        logging.logging.info("building fortran binaries")
        command = ["make", "ddscat"]
        if subprocess.call(command) != 0:
            sys.exit(-1)
        if self.bdist_dir:
            logging.logging.info("Copying executables...")
            output_dir = self.bdist_dir / self.pkg_name
            output_dir.mkdir(parents=True, exist_ok=True)
            for executable_name in [
                "ddscat",
                "ddpostprocess",
                "vtrconvert",
                "calltarget",
            ]:
                shutil.copy2(Path(f"ddscat/src/{executable_name}"), output_dir)


class CustomBuild(build):
    sub_commands = [("build_custom", None)] + build.sub_commands


setup(cmdclass={"build": CustomBuild, "build_custom": CustomCommand})
