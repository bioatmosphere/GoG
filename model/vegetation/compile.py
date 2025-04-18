import os
import subprocess

import numpy as np


os.chdir("src/builddir/")

# First, compile the Fortran code with f2py
constants_f90 = "Constants.f90"
compile_command = f"python -m numpy.f2py -c {constants_f90} -m constants --backend meson"



print("Compiling Fortran module...")
result = subprocess.run(compile_command, shell=True, capture_output=True, text=True)
if result.returncode != 0:
    print("Compilation failed:")
    print(result.stderr)
    exit(1)
print("Compilation successful!")