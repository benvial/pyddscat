#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/pyddscat
"""
Functions for reading and writing from files.

"""

import subprocess
import os
import os.path
import time
import numpy as np
import posixpath
import pdb
import copy

from . import core
from . import utils
from . import ranges
from . import targets
from .__about__ import __version__


def build_ddscat_par(settings, target):
    """
    Return a string with the contents of the ddscat.par file.

    :param settings: a :class:`core.Settings` object
    :param target: a :class:`core.Target` object
    """
    out = ""
    out += f"=== Generated by pyddscat v{__version__} ({time.asctime()}) ===\n"

    out += "**** Preliminaries ****\n"
    out += "DOTORQ\n" if settings.CMDTRQ else "NOTORQ\n"
    out += settings.CMDSOL + "\n"
    out += settings.CMDFFT + "\n"
    out += settings.CALPHA + "\n"
    out += settings.CBINFLAG + "\n"

    out += "**** Initial Memory Allocation ****" + "\n"
    if settings.InitialMalloc is not None:
        out += settings.InitialMalloc.__str__()[1:-1] + "\n"
    elif isinstance(target, targets.FROM_FILE):
        out += str(target.sh_param())[1:-1] + "\n"
    else:
        out += "100 100 100\n"

    out += target.save_str()  # Target defn goes here

    out += "**** Additional Nearfield calculation? ****\n"
    out += "1\n" if settings.NRFLD else "0\n"
    out += settings.NRFLD_EXT.__str__()[1:-1] + "\n"

    out += "**** Error Tolerance ****\n"
    out += str(settings.TOL) + "\n"

    out += "**** Maximum number of iterations allowed ****\n"
    out += str(settings.MXITER) + "\n"

    out += "**** Interaction cutoff parameter for PBC calculations ****\n"
    out += str(settings.GAMMA) + "\n"

    out += "**** Angular resolution for calculation of <cos>, etc. ****\n"
    out += str(settings.ETASCA) + "\n"

    out += "**** Vacuum wavelengths (micron) ****\n"
    out += settings.wavelengths.__str__() + "\n"

    out += "**** Refractive index of ambient medium\n"
    out += str(settings.NAMBIENT) + "\n"

    out += "**** Effective Radii (micron) **** \n"
    if isinstance(settings.scale_range, ranges.How_Range):
        scale = copy.copy(settings.scale_range)
    else:
        scale = ranges.How_Range(settings.scale_range, settings.scale_range, 1)

    aeff = scale
    aeff.first *= target.aeff
    aeff.last *= target.aeff
    out += str(aeff) + "\n"

    out += "**** Define Incident Polarizations ****\n"
    out += utils.complexV2str(settings.Epol) + "\n"
    out += "2\n" if settings.IORTH else "1\n"

    out += "**** Specify which output files to write ****\n"
    out += "1\n" if settings.IWRKSC else "0\n"

    out += "**** Prescribe Target Rotations ****\n"
    out += settings.beta.__str__() + "\n"
    out += settings.theta.__str__() + "\n"
    out += settings.phi.__str__() + "\n"

    out += "**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****\n"
    out += settings.initial.__str__()[1:-1] + "\n"

    out += "**** Select Elements of S_ij Matrix to Print ****" + "\n"
    out += str(len(settings.S_INDICES)) + "\n"
    for s in settings.S_INDICES:
        out += "%d " % s
    out += "\n"

    out += "**** Specify Scattered Directions ****\n"
    out += settings.CMDFRM + "\n"
    out += str(len(settings.scat_planes)) + "\n"
    for s in settings.scat_planes:
        out += str(s) + "\n"
    out += "\n"

    return out


def _parseline(line):
    """
    Process a line from the DDSCAT file.

    :param line: The input string to process
    :returns: A string with extraneous characters removed

    Ignores any characters after a '=' or '!'
    Removes quote characters
    """

    if pts := [line.find(c) for c in "=!" if line.find(c) != -1]:
        line = line[: min(pts)]

    # Remove ' and "
    # line = line.translate(None, '\'\"')
    line = line.translate(str.maketrans("", "", "'\""))

    # Remove leading and trailing whitespace
    line = line.strip()

    return line


def QSub_Batchfile(fname, base_path, folders):
    """
    Create a csh script for batch submission of many runs via qsub.

    This assumes that the server uses posix paths, regardless of the path
    convention on the local machine.

    :param fname: the name of the batch file
    :param base_path: the path from which the folders will be resolved.
        This must be an absolute path on the server.
    :param folders: a list of folders (relative to base_path) containing
                 the submission scripts (.sge files)
    """

    norm = posixpath.normpath
    join = posixpath.join

    with open(fname, "wt") as f:
        f.write("#!/bin/csh\n")
        for l in folders:
            folder = norm(join(base_path, norm(l)))
            sge_file = join(folder, "submit.sge")
            f.write("qsub -wd %s %s \n" % (folder, sge_file))

    try:
        subprocess.call(["chmod", "+x", fname])
    except OSError:
        pass


def build_constant_material(index, fname):
    index = complex(index)
    re, im = index.real, index.imag
    contents = f"""Custom Constant Material
    1 2 3 0 0 = columns for wave, Re(n), Im(n)
    wave(um) Re(n)  Im(n)
    0.0   {re}  {im}
    1.0   {re}  {im}
    2.0   {re}  {im}
    """
    with open(fname, "wt") as f:
        f.write(contents)


class Line:
    """Class representing a line."""

    def __init__(self, start, end, N):
        self.start = start
        self.end = end
        self.N = N


def generate_ddpostprocess_par(
    infile="w000r000k000.E1", prefix="VTRoutput", vtr=False, lines=[]
):
    line = lines != []
    tpl = f"""'{infile}'
'{prefix}'
{int(vtr)}
{int(line)}"""
    for line in lines:
        s, e = line.start, line.end
        tpl += f"\n{s[0]} {s[1]} {s[2]} {e[0]} {e[1]} {e[2]} {line.N}"

    return tpl


def write_ddpostprocess_par(fname, *args, **kwargs):
    tpl = generate_ddpostprocess_par(*args, **kwargs)
    with open(fname, "wt") as f:
        f.write(tpl)
