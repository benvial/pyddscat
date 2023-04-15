#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/pyddscat
"""
A set of tools for setting up and working with the program DDScat.
"""

import tempfile
import subprocess
import warnings
import numpy as np
import os
import os.path
import posixpath, ntpath
import copy
import pdb
import inspect

from . import targets
from . import results
from . import fileio
from . import ranges
from . import utils

#: The configuration settings
config = {}

#: Define polarization states using the spectroscopists' convention
pol_cL = np.array([0, 1 + 0j, 0 + 1j])
pol_cR = np.array([0, 0 + 1j, 1 + 0j])
pol_lH = np.array([0, 1 + 0j, 0 + 0j])
pol_lV = np.array([0, 0 + 0j, 1 + 0j])

##: Define polarization states using DDSCAT's convention
# pol_cR=np.array([0, 1+0j, 0+1j])
# pol_cL=np.array([0, 0+1j, 1+0j])
# pol_lH=np.array([0, 1+0j, 0+0j])
# pol_lV=np.array([0, 0+0j, 1+0j])


class Settings(object):
    """
    DDSCAT execution parameters

    Most of the field names correspond to their definitions in the ddscat.par file.
    Details of the target are stored in the ```Target``` definition, not here.

    """

    #: Either do or skip torque calculations (True:'DOTORQ', False:'NOTORQ')
    CMDTRQ = False

    #: Solution method (PBCGS2, PBCGST, GPBICG, PETRKP, QMRCCG)
    CMDSOL = "PBCGS2"  #: = CMDSOL*6 (PBCGS2, PBCGST, GPBICG, PETRKP, QMRCCG) --

    #: FFT method (GPFAFT, FFTMKL)
    CMDFFT = "GPFAFT"

    #: Prescription for polarizabilities (GKDLDR, LATTDR)
    CALPHA = "GKDLDR"

    #: Specify binary output
    CBINFLAG = "NOTBIN"

    #: Initial Memory Allocation. None will attempt to automatically size the initial malloc
    InitialMalloc = None

    #: Either do or skip nearfield calculations (True, False)
    NRFLD = False

    #: Fractional extension of calculated volume in (-x,+x,-y,+y,-z,+z)
    NRFLD_EXT = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    #: Error Tolerence.
    TOL = 1.00e-5

    #: Maximum number of iterations
    MXITER = 600

    #: Interaction cutoff parameter for PBC calculations (1e-2 is normal, 3e-3 for greater accuracy)
    GAMMA = 1.00e-2

    #: Angular resolution for calculation of <cos>, etc. (number of angles is proportional to [(3+x)/ETASCA]^2 )
    ETASCA = 0.5

    #: Specify which output files to write (True: write ".sca" file for each target orient. False: suppress)
    IWRKSC = True

    #: Vacuum wavelengths (first,last,how many,how=LIN,INV,LOG)
    wavelengths = ranges.How_Range(0.3500, 0.8000, 10, "LIN")

    #: Refractive index of ambient medium
    NAMBIENT = 1.000

    #: Define a range of scales for the particle geometry. A number indicates a single size calc.
    scale_range = ranges.How_Range(1, 1, 1)

    #: Define Incident Polarizations (Polarization state e01 (k along x axis)
    Epol = pol_lV

    #: Specify whether to calculate orthogonal polarization state (True, False)
    IORTH = True

    #: Prescribe Target Rotation beta (rotation around a1). (betamin, betamx, nbeta)
    beta = ranges.Lin_Range(0.0, 0.0, 1)

    #: Prescribe Target Rotation theta (angle between a1 and k). (thetamin, tetamx, ntheta)
    theta = ranges.Lin_Range(0.0, 0.0, 1)

    #: Prescribe Target Rotation phi (rotation of a1 around k). (phimin, phimax, nphi)
    phi = ranges.Lin_Range(0.0, 0.0, 1)

    #: Specify first IWAV, IRAD, IORI (0 0 0 to begin fresh)
    initial = [0, 0, 0]

    #: Select Elements of S_ij Matrix to Print
    S_INDICES = [11, 12, 13, 14, 21, 22, 31, 41, 44]

    #: Specify reference frame for scattering ('LFRAME', 'TFRAME')
    CMDFRM = "LFRAME"

    #: Specify Scattered Directions
    scat_planes = [ranges.Scat_Range(0, 0, 180, 5), ranges.Scat_Range(90, 0, 180, 5)]

    def __init__(self, folder=None, **kwargs):
        # If available, settings come from default.par file
        default = utils.resolve_profile("default.par")
        if default is not None:
            kwargs = dict(
                list(self._read_values(default).items()) + list(kwargs.items())
            )

        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(f"Settings has no attribute {key}")

    @classmethod
    def fromfile(cls, fname):
        """
        Create new Settings object with values loaded from the specified file.
        """

        return cls(**cls._read_values(fname))

    @classmethod
    def _read_values(cls, fname):
        """
        Load the values for a new Settings object from the specified file.

        :param fname: The filename of the ddscat.par file to loadfrom
        :returns: A dict of the values
        """

        with open(fname, "r") as f:
            lines = [fileio._parseline(line) for line in f.readlines()]
        # Process target.
        # This is discarded, only needed to determine line spacing and
        # style of scattering specifier
        n_mat = int(lines[12])
        target = targets.Target.fromfile(fname)
        del lines[9 : 13 + n_mat]

        # Process settings
        settings = {
            "CMDTRQ": lines[2] == "DOTORQ",
            "CMDSOL": lines[3],
            "CMDFFT": lines[4],
            "CALPHA": lines[5],
            "CBINFLAG": lines[6],
            "InitialMalloc": np.array(list(map(int, lines[8].split()))),
        }
        settings["NRFLD"] = bool(int(lines[10]))
        settings["NRFLD_EXT"] = np.array(list(map(float, lines[11].split()[:6])))

        settings["TOL"] = float(lines[13])
        settings["MXITER"] = int(lines[15])
        settings["GAMMA"] = float(lines[17])
        settings["ETASCA"] = float(lines[19])
        settings["wavelengths"] = ranges.How_Range.fromstring(lines[21])

        settings["NAMBIENT"] = float(lines[23])

        scale = ranges.How_Range.fromstring(lines[25])
        scale.last /= scale.first
        scale.first = 1.0
        settings["scale_range"] = scale

        settings["Epol"] = utils.str2complexV(lines[27])
        settings["IORTH"] = int(lines[28]) == 2

        settings["IWRKSC"] = bool(int(lines[30]))

        settings["beta"] = ranges.Lin_Range.fromstring(lines[32])
        settings["theta"] = ranges.Lin_Range.fromstring(lines[33])
        settings["phi"] = ranges.Lin_Range.fromstring(lines[34])

        if lines[36].find(",") != -1:
            settings["initial"] = list(map(int, lines[36].split(",")))
        else:
            settings["initial"] = list(map(int, lines[36].split()))

        settings["S_INDICES"] = list(map(int, lines[39].split()))

        settings["CMDFRM"] = lines[41]
        n_scat = int(lines[42])

        if n_scat > 0:
            if isinstance(target, targets.Periodic):
                if target.dimension == 1:
                    settings["scat_planes"] = [
                        ranges.Scat_Range_1dPBC.fromstring(line)
                        for line in lines[43 : 43 + n_scat]
                    ]
                elif target.dimension == 2:
                    settings["scat_planes"] = [
                        ranges.Scat_Range_2dPBC.fromstring(line)
                        for line in lines[43 : 43 + n_scat]
                    ]

            else:
                settings["scat_planes"] = [
                    ranges.Scat_Range.fromstring(line)
                    for line in lines[43 : 43 + n_scat]
                ]
        return settings

    def copy(self):
        """
        Make a copy of these settings.
        """
        return copy.deepcopy(self)


class DDscat(object):
    """
    A class for managing a DDSCAT run.

    Loosely a DDscat object includes two parts: a settings file and a target.

    All parameters are optional. If they are not specified, DDscat defaults
    to the settings found in the file default.par.
    :param folder: The subfolder in which to store files.
    :param settings: The :class:`Settings` object for this run.
    :param target: The :class:`targets.Target` for this run. Defaults to a Au sphere.


    Example::

        # Build a target
        t=CYLNDRCAP(0.100, 0.030):

        # Initialize a DDscat run
        d=DDscat(target=t)

        # Modify run parameters as desired
        d.settings.NAMBIEND=1.5
        d.settings.Epol=np.array([0, 1j, 1])

        # Write the output and run the simulation
        d.calculate()

        # Load the results
        r=results.FolderCollection()

        # Plot them
        r.plot()
    """

    def __init__(self, folder=None, settings=None, target=None):
        if folder is None:
            self._folder = tempfile.mkdtemp()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
            self._folder = folder

        self.settings = Settings() if settings is None else settings.copy()
        self.settings.folder = self._folder
        if target is None:
            self.target = targets.Target()
        else:
            self._target = target

        self.target.folder = self._folder

    @classmethod
    def fromfile(cls, fname):
        """
        Create new DDscat object with values loaded from the specified file.
        """
        settings = Settings.fromfile(fname)
        target = targets.Target.fromfile(fname)
        return cls(settings=settings, target=target)

    def __str__(self):
        """
        A string representation.
        This should be expanded
        """
        return "DDScat Definition"

    def copy(self):
        """
        Make a copy of this run.
        """
        return copy.deepcopy(self)

    def write(self, *args, **kwargs):
        """Write the .par file and target definitions to file.

        :param args: Optional arguments and kwargs are passed to the profile_script

        :param submit_script: Use True to write a .sge file for submitting the job
                            to an SGE cluster via ```qsub```.

        If the current profile defines a function ```write_script``` then this
        is run after the ddscat.par has been written, and is called with any
        addition arguments to write.
        """

        self.check()

        s = fileio.build_ddscat_par(self.settings, self.target)

        with open(os.path.join(self.folder, "ddscat.par"), "wt") as f:
            f.write(s)

        self.target.folder = self.folder
        self.target.write()

        try:
            config["write_script"](self, config, *args, **kwargs)
        except KeyError:
            pass

    def check(self):
        """
        Perform simple self consistency checks on the job.

        Raises an exception if it finds a problem.
        """

        # Check that the correct type of scattering plane is used for Periodic targets
        if isinstance(self.target, targets.Periodic):
            for s in self.settings.scat_planes:
                if self.target.dimension == 1:
                    if not isinstance(s, ranges.Scat_Range_1dPBC):
                        raise TypeError(
                            "Target is 1D periodic. Scattering plane definitions must be Scat_Range_1dPBC"
                        )
                elif self.target.dimension == 2:
                    if not isinstance(s, ranges.Scat_Range_2dPBC):
                        raise TypeError(
                            "Target is 2D periodic. Scattering plane definitions must be Scat_Range_2dPBC"
                        )

            if self.settings.CMDFRM == "LFRAME":
                raise ValueError("For PBC calculations CMDFRM must be 'TFRAME'")

        else:
            for s in self.settings.scat_planes:
                if not isinstance(s, ranges.Scat_Range):
                    raise TypeError(
                        "Target is finite, isolated. Scattering plane definitions must be Scat_Range"
                    )

    @property
    def folder(self):
        """This run's home folder"""
        return self._folder

    @folder.setter
    def folder(self, newfolder):
        """Redefine the run's home folder."""
        if not os.path.exists(newfolder):
            os.makedirs(newfolder)

        self._folder = newfolder
        self.target.folder = newfolder

    @property
    def target(self):
        """The run's target"""
        return self._target

    @target.setter
    def target(self, newtarget):
        """Redefine the run's target."""
        self._target = newtarget
        self._target.folder = self.folder

    def info(self):
        """Print some basic run info"""
        wave = self.settings.wavelengths.table

        table = np.vstack([wave, self.mkd, self.x, self.alpha, self.beta])
        print("Target: ", self.target.directive, "(", self.target.__class__, ")")
        print("N=%d" % self.target.N)
        print("aeff=%f" % self.target.aeff)
        print("d=%f" % self.target.d)
        print("wave          mkd         x           alpha       beta")
        print("--------------------------------------------------------------")
        for r in table.transpose():
            print(r)

    @property
    def x(self):
        """Calculate the x-parameter (Userguide p8)."""
        a = self.target.aeff

        out = [2 * np.pi * a / wl for wl in self.settings.wavelengths]

        return np.asarray(out)

    @property
    def mkd(self):
        """Calculate m*k*d (Userguide p8)."""
        m_dat = results.MInTable(self.target.material[0])

        out = []
        k = 2 * np.pi * self.target.d
        for wl in self.settings.wavelengths:
            m = m_dat(wl)
            m = np.abs(m["Rem"] + 1j * m["Imm"])
            out.append(k / wl * m)

        return np.asarray(out)

    @property
    def alpha(self):
        """Calculate the alpha parameter (Userguide Eqn 8, p8)."""
        N = self.target.N
        m_dat = results.MInTable(self.target.material[0])

        out = []
        for wl in self.settings.wavelengths:
            m = m_dat(wl)
            k = 9.88 * wl * (N / 10**6) ** (1 / 3)
            out.append(k / np.abs(m["Rem"] + 1j * m["Imm"]))

        return np.asarray(out)

    @property
    def beta(self):
        """Calculate the beta parameter (Userguide Eqn 8, p8)."""
        N = self.target.N
        m_dat = results.MInTable(self.target.material[0])
        out = []
        for wl in self.settings.wavelengths:
            m = m_dat(wl)
            out.append(62 / np.abs(m["Rem"] + 1j * m["Imm"]) * (N / 10**6) ** (1 / 3))

        return np.asarray(out)

    def calculate(self, silent=True):
        """Start local calculation of this run.

        This assumes that ```ddscat``` is in the path.

        :param silent: If true suppresses output to the screen
        """

        self.write()

        command = os.path.join(config["ddscat_path"], "ddscat")

        # try:
        #     __IPYTHON__
        # except NameError:
        #     pass
        # else:
        #     silent=True
        #     warnings.warn('DDscat output will not be displayed in IPython')

        if silent:
            subprocess.call(f"{command} 2> output.log", shell=True, cwd=self.folder)
        else:
            utils.logger.debug("Starting calculation...")
            subprocess.call(
                f"{command} 2>&1 | tee output.log", shell=True, cwd=self.folder
            )
            utils.logger.debug("Done!")

    def VTRconvert(self, outfile=None):
        """
        Convert the target shape into a form viewable in Paraview

        """
        self.target.VTRconvert()

    def calltarget(self):
        """
        Executes ```calltarget``` to generate a target.out file for builtin target geometries.
        """

        self.write()

        subprocess.call(
            os.path.join(config["ddscat_path"], "calltarget"), cwd=self.folder
        )

    def postprocess_fields(
        self,
        infile="w000r000k000.E1",
        prefix="VTRoutput",
        vtr=False,
        lines=[],
        verbose=False,
    ):
        """Execute ddpostprocess"""

        fname = os.path.join(self.folder, "ddpostprocess.par")
        fileio.write_ddpostprocess_par(fname, infile, prefix, vtr, lines)
        result = subprocess.run(
            ["ddpostprocess"],
            cwd=self.folder,
            capture_output=True,
            text=True,
        )
        if verbose:
            print(result.stdout)
            print(result.stderr)
        datafile = os.path.join(self.folder, "ddpostprocess.out")
        data = np.loadtxt(datafile, skiprows=1)
        self.fields = dict(
            x=data[:, 0],
            y=data[:, 1],
            z=data[:, 2],
            Ex=data[:, 3] + 1j * data[:, 4],
            Ey=data[:, 5] + 1j * data[:, 6],
            Ez=data[:, 7] + 1j * data[:, 8],
        )
        return self.fields

    @property
    def output_collection(self):
        return results.FolderCollection(path=self.folder)

    @property
    def output(self):
        return self.output_collection[self.folder]

    # @property
    # def fields(self):
    #     fname = os.path.join(self.folder, "w000r000k000.E1")
    #     return results.EnTable(fname)


def set_config(fname=None):
    """
    Select which configuration profile to use for pyddscat

    :param fname: The name of the file that contains the configuration
                  If None then load the default profile

    Profiles are stored in Python script files.

    The search scheme is to first look for the file in the CWD, followed by
    the folder ~/.pyddscat/ and finally the subdiretory profiles/ relative to
    where the config.py module resides.
    """
    global config

    if fname is None:
        fname = "default.py"

    if not fname.endswith(".py"):
        fname += ".py"

    full_name = utils.resolve_profile(fname)

    if full_name is None:
        raise IOError

    exec(compile(open(full_name).read(), full_name, "exec"), {}, config)

    # Remove any imported modules from config
    for k, v in list(config.items()):
        if inspect.ismodule(v):
            del config[k]

    config["profile"] = os.path.abspath(full_name)

    # Associate the correct path style based on OS
    if config["os"].lower() in ["unix", "mac"]:
        config["path_style"] = posixpath
    elif config["os"].lower() == "windows":
        config["path_style"] = ntpath
    else:
        raise ValueError(f"Unknown OS: {config['os']}")


# ON import set the configuration to default
set_config(None)
