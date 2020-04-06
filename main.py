import excitons
import reader_cif
import reader_orca
import reader_cube
import numpy as np
import linker
from os import getcwd
import kujo_io
import concurrent.futures
from threading import active_count


def integrate_translated_cube(vars: list):
    is_bohr = None
    try:
        translate = np.array([float(vars[0]), float(vars[1]), float(vars[2])])
        if len(vars) > 3:
            is_bohr = bool(vars[3])
    except ValueError:
        kujo_io.output_error("Value type error at options!", -1)
    if is_bohr is None:




dispatcher = {
    "integrate_translated_cube": integrate_translated_cube
}

threads_number = active_count()
cube = reader_cube.Cube("placeholder")
instructions, options = kujo_io.read_input("input.txt")
for i in range(len(instructions)):
    dispatcher[instructions[i]](options)