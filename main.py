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


def options_parse(dsp: dict, opt: list):
    for i1 in range(0, len(opt), 2):
        if opt[i1] in dsp:
            if opt[i1 + 1].lower() == "true":
                dsp[opt[i1]] = True
            elif opt[i1 + 1].lower() == "false":
                dsp[opt[i1]] = False
            else:
                try:
                    fl = float(opt[i1 + 1])
                    dsp[opt[i1]] = fl
                except ValueError:
                    if "[" in opt[i1 + 1] or "]" in opt[i1 + 1]:
                        opt[i1 + 1].replace("[", "")
                        opt[i1 + 1].replace("]", "")
                        words = opt[i1 + 1].split(",")
                        temp = np.zeros([len(words), 1])
                        try:
                            for i in range(len(words)):
                                temp[i, 0] = float(words[i])
                            dsp[opt[i1]] = temp
                        except ValueError:
                            exit(-1)
                    else:
                        dsp[opt[i1]] = opt[i1 + 1]
        else:
            exit(-1)


def read_cube(v: list):
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cube
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cube = reader_cube.Cube(full_path)


def read_cif(v: list):
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cif
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cif = reader_cif.CifFile(full_path)


def integrate_translated_cube(v: list):
    options_dispatcher = {
        "is_bohr": False,
        "vector": None,
        "vector_cif": "",
        "multiplier": 1.0,
    }
    options_parse(options_dispatcher, v)
    low_limits = []
    up_limits = []
    low_limits.append(np.array([0, cube.steps[1, 0], cube.steps[2, 0]]))
    step = np.array([cube.steps[0, 0] // threads_number, 0, 0])
    for _ in range(threads_number - 1):
        up_limits.append(low_limits[-1] + step)
        low_limits.append(up_limits[-1])
    up_limits.append(np.array([cube.steps[0, 0], cube.steps[1, 0], cube.steps[2, 0]]))
    if options_dispatcher["vector_cif"] == "a":
        translate = cif.vector_a * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "b":
        translate = cif.vector_b * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "c":
        translate = cif.vector_c * options_dispatcher["multiplier"]
    elif options_dispatcher["vector"] is not None:
        translate = options_dispatcher["vector"] * options_dispatcher["multiplier"]
    else:
        exit(-1)
    if options_dispatcher["is_bohr"] is False:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(cube.integrate, low_limits[i], up_limits[i], translate)
                       for i in range(threads_number)]
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(cube.integrate, low_limits[i], up_limits[i], translate, True)
                       for i in range(threads_number)]
    r = 0
    for x in results:
        r = r + x.result()
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
        result = dispatcher[instructions[i]](options[i])
        full_path = getcwd() + "/output.txt"
        if result is not None:
            file = open(full_path, "a+")
            file.write(f"Instruction: '{instructions[i]}' at the line {i} returned: {result}\n")
            file.close()
