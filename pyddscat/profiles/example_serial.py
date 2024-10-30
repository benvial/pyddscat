#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat

# The profile name
name = "luna_serial"

# The operating system
os = "unix"

# The absolute path to the folder containing materials properties
mat_library = "/home/guests/mark/mat_prop"

# The absolute path of the ddscat executable
mat_library = "/cluster/bin/ddscat"


def write_script(job, config):
    import os.path

    with open(os.path.join(job.folder, "submit.sge"), "wt") as f:
        f.write("#!/bin/csh\n")
        f.write("#\n#\n#\n")
        f.write("# ---------------------------\n")
        f.write("# our name \n")

        f.write("#$ -N ddscat_ser_PRE_\n")
        f.write("#\n")

        f.write("# stderr >& stdout\n")
        f.write("#$ -j y\n")
        f.write("#\n")
        f.write("# ---------------------------\n")

        f.write("set hostname=`/bin/hostname`\n")

        f.write("echo beginning `pwd`\n")
        f.write("date\n")
        f.write("time /cluster/bin/ddscat\n" % config["ddscat_exec"])
        f.write("echo completed `pwd`\n")
        f.write("echo '------------------------------------------'\n")
        f.write("date\n")

        f.write(
            'foreach old (ddscat_ser_PRE_.*)\nmv $old "$old:gas/_PRE_//"".txt"\nend\n'
        )
