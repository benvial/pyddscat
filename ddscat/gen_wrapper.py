import sys

execname = sys.argv[1]


code = f"""#!/usr/bin/env python

import subprocess, sys
import os.path

cmd = os.path.join(os.path.dirname(__file__), "{execname}_bin")
subprocess.call([cmd] + sys.argv[1:])
"""


with open(f"{execname}", "w") as f:
    f.write(code)
