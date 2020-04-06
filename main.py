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
        pass
    low_limits = []
    up_limits = []
    low_limits.append(np.array([0, 0, 0]))
    step = np.array([cube.steps[0, 0] // threads_number, cube.steps[1, 0], cube.steps[2, 0]])
    for _ in range(threads_number - 1):
        up_limits.append(low_limits[-1] + step)
        low_limits.append(up_limits[-1])
    up_limits.append(np.array([cube.steps[0, 0], cube.steps[1, 0], cube.steps[2, 0]]))
    if is_bohr is None:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(cube.integrate, low_limits[i], up_limits[i], translate)
                       for i in range(threads_number)]
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(cube.integrate, low_limits[i], up_limits[i], translate, is_bohr)
                       for i in range(threads_number)]
    r = 0
    for x in results:
        r = r + x
    return r



dispatcher = {
    "integrate_translated_cube": integrate_translated_cube
}

threads_number = active_count()
cube = reader_cube.Cube("placeholder")
instructions, options = kujo_io.read_input("input.txt")
for i in range(len(instructions)):
    dispatcher[instructions[i]](options)