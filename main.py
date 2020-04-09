import numpy as np
from os import getcwd
import sys
import kujo_io
import concurrent.futures


cube = None
cif = None
cluster = None
hard_cutoff = 0.0
int_cutoff = 0.0


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
    if "reader_cube" not in sys.modules:
        import reader_cube
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cube
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cube = reader_cube.Cube(full_path)


def read_cif(v: list):
    if "reader_cif" not in sys.modules:
        import reader_cif
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cif
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cif = reader_cif.CifFile(full_path)


def integrate_translated_cube(v: list):
    if "reader_cube" not in sys.modules:
        import reader_cube
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


def set_hard_cutoff(v: list):
    options_dispatcher = {
        "value": ""
    }
    options_parse(options_dispatcher, v)
    global hard_cutoff
    hard_cutoff = float(options_dispatcher["value"])


def set_int_cutoff(v: list):
    options_dispatcher = {
        "value": ""
    }
    options_parse(options_dispatcher, v)
    global int_cutoff
    int_cutoff = float(options_dispatcher["value"])


def build_cluster(v: list):
    if "reader_cif" not in sys.modules:
        import reader_cif
    global cluster
    options_dispatcher = {
        "a": 0,
        "b": 0,
        "c": 0
    }
    options_parse(options_dispatcher, v)
    cluster = reader_cif.Cluster(cif)
    cluster.build()
    it = list(range((len(cluster.pre_molecules) * cluster.pre_molecules[0].num_atoms)))
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(cluster.connectivity, it)
    cluster.rebuild()
    it = []
    for i in range(len(cluster.pre_molecules)):
        it.append(cluster.pre_molecules[i])
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(cluster.to_cartesian, it)
    cluster.multiply(options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])


dispatcher = {
    "integrate_translated_cube": integrate_translated_cube,
    "read_cube": read_cube,
    "read_cif": read_cif,
    "build_cluster": build_cluster,
    "print_cluster_xyz": cluster.print_cluster_xyz,
    "print_cluster_kujo": cluster.print_cluster_readable,
    "set_hard_cutoff": set_hard_cutoff,
    "set_int_cutoff": set_int_cutoff
}

if __name__ == "__main__":
    instructions, options = kujo_io.read_input("input.txt")
    out = "/output.out"
    for i in range(len(instructions)):
        if len(options) != 0:
            result = dispatcher[instructions[i]](options[i])
        else:
            result = dispatcher[instructions[i]]()
        full_path = getcwd() + out
        if result is not None:
            file = open(full_path, "a+")
            file.write(f"{instructions[i]} with options {options[i]} returned: {result}\n")
            file.close()
