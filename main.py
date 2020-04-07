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


cube = None
cif = None


def read_cube(v: list):
    global cube
    full_path = getcwd() + f"/{v[0]}"
    cube = reader_cube.Cube(full_path)


def read_cif(v: list):
    global cif
    full_path = getcwd() + f"/{v[0]}"
    cif = reader_cif.CifFile(full_path)


def integrate_translated_cube(v: list):
    is_bohr = None
    try:
        translate = np.array([float(v[0]), float(v[1]), float(v[2])])
        if len(v) > 3:
            is_bohr = bool(v[3])
    except ValueError:
        kujo_io.output_error("Type error at options!", -1)
        pass
    low_limits = []
    up_limits = []
    low_limits.append(np.array([0, 0, 0]))
    step = np.array([сube.steps[0, 0] // threads_number, cube.steps[1, 0], cube.steps[2, 0]])
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
    "integrate_translated_cube": integrate_translated_cube,
    "read_cube": read_cube,
    "read_cif": read_cif
}

if __name__ == "__main__":
    threads_number = active_count()
    instructions, options = kujo_io.read_input("input.txt")
    for i in range(len(instructions)):
        result = dispatcher[instructions[i]](options)
        full_path = getcwd() + "/output.txt"
        if result is not None:
            file = open(full_path, "a+")
            file.write(f"Instruction: '{instructions[i]}' at the line {i} returned: {result}\n")
            file.close()
