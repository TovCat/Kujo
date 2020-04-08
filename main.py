#import excitons
import reader_cif
import reader_orca
import reader_cube
import numpy as np
import linker
from os import getcwd
import kujo_io
import concurrent.futures


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
        "vector": None,
        "vector_cif": "",
        "multiplier": 1.0,
    }
    options_parse(options_dispatcher, v)
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
    r = 0.0
    it = []
    for i1 in range(cube.steps[0, 0]):
        it.append([translate, i1])
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(cube.integrate, it)
        for x in results:
            r = r + x
    return r


dispatcher = {
    "integrate_translated_cube": integrate_translated_cube,
    "read_cube": read_cube,
    "read_cif": read_cif
}

if __name__ == "__main__":
    instructions, options = kujo_io.read_input("input.txt")
    for i in range(len(instructions)):
        result = dispatcher[instructions[i]](options[i])
        full_path = getcwd() + "/output.txt"
        if result is not None:
            file = open(full_path, "a+")
            file.write(f"Instruction: '{instructions[i]}' at the line {i} returned: {result}\n")
            file.close()
