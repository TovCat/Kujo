import numpy as np
from os import getcwd
from copy import deepcopy
import utility.kujo_io
import concurrent.futures
import reader.cube
import reader.charges
import reader.orca
import reader.cif
import time
import energy.excitons

cube = None
cif = None
cluster = None
charges = None
orca = None
H = None
disorders = []
max_w = None
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
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cube
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cube = reader.cube.Cube(full_path)


def read_cif(v: list):
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global cif
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cif = reader.cif.CifFile(full_path)


def read_charges(v: list):
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global charges
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    charges = reader.charges.Charges(full_path)


def read_orca(v: list):
    options_dispatcher = {
        "file": ""
    }
    options_parse(options_dispatcher, v)
    global orca
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    orca = reader.orca.Orca(full_path)


def translated_coupling_td_integration(v: list):
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
    if max_w is not None:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_w) as executor:
            results = executor.map(cube.integrate, it)
            for x in results:
                r = r + x
    else:
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


def set_max_workers(v: list):
    options_dispatcher = {
        "value": ""
    }
    options_parse(options_dispatcher, v)
    global max_w
    max_w = int(options_dispatcher["value"])


def build_cluster(v: list):
    global cluster
    global cif
    options_dispatcher = {
        "a": 0,
        "b": 0,
        "c": 0
    }
    options_parse(options_dispatcher, v)
    cluster = reader.cif.Cluster(cif, options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])
    cluster.build()
    it = []
    cluster.to_cartesian()
    cluster.connectivity()
    cluster.rebuild()
    cluster.find_rotations()
    cluster.multiply(options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])
    cluster.build_rmc()


def translated_coupling_extended_dipole(v: list):
    options_dispatcher = {
        "vector": None,
        "vector_cif": "",
        "multiplier": 1.0,
        "d": 0.0,
    }
    options_parse(options_dispatcher, v)
    if options_dispatcher["vector_cif"] == "a":
        t = cif.vector_a * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "b":
        t = cif.vector_b * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "c":
        t = cif.vector_c * options_dispatcher["multiplier"]
    elif options_dispatcher["vector"] is not None:
        t = options_dispatcher["vector"] * options_dispatcher["multiplier"]
    else:
        exit(-1)
    mu = orca.mu
    mol1 = cluster.molecules[0]
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    d = options_dispatcher["d"]
    return energy.excitons.coupling_extended_dipole(mol1, mol2, mu, d)


def translated_coupling_dipole(v: list):
    options_dispatcher = {
        "vector": None,
        "vector_cif": "",
        "multiplier": 1.0,
    }
    options_parse(options_dispatcher, v)
    global cif
    if options_dispatcher["vector_cif"] == "a":
        t = cif.vector_a * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "b":
        t = cif.vector_b * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "c":
        t = cif.vector_c * options_dispatcher["multiplier"]
    elif options_dispatcher["vector"] is not None:
        t = options_dispatcher["vector"] * options_dispatcher["multiplier"]
    else:
        exit(-1)
    mu = orca.mu
    mol1 = cluster.molecules[0]
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    r = mol1.mass_center() - mol2.mass_center()
    return energy.excitons.coupling_dipole(mol1, mol2, mu)


def translated_coupling_charges(v: list):
    options_dispatcher = {
        "vector": None,
        "vector_cif": "",
        "multiplier": 1.0,
    }
    options_parse(options_dispatcher, v)
    if options_dispatcher["vector_cif"] == "a":
        t = cif.vector_a * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "b":
        t = cif.vector_b * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "c":
        t = cif.vector_c * options_dispatcher["multiplier"]
    elif options_dispatcher["vector"] is not None:
        t = options_dispatcher["vector"] * options_dispatcher["multiplier"]
    else:
        exit(-1)
    mol1 = charges.mol
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    return energy.excitons.coupling_charges(mol1, mol2, charges.q)


def generate_disorder(v: list):
    options_dispatcher = {
        "sigma": 0,
        "n": 1
    }
    options_parse(options_dispatcher, v)
    global disorders
    n = H.shape[0]
    for i in range(options_dispatcher["n"]):
        disorders.append(n, options_dispatcher["sigma"])


def calculate_coupling(v: list):
    options_dispatcher = {
        "method": "",
        "site1": 0,
        "site2": 0
    }
    options_parse(options_dispatcher, v)
    mol1 = cluster.molecules[options_dispatcher["site1"]]
    mol2 = cluster.molecules[options_dispatcher["site2"]]
    if options_dispatcher["method"] == "dipole":
        return
    

dispatcher = {
    "translated_coupling_td_integration": translated_coupling_td_integration,
    "translated_coupling_extended_dipole": translated_coupling_extended_dipole,
    "translated_coupling_dipole": translated_coupling_dipole,
    "translated_coupling_charges": translated_coupling_charges,
    "read_cube": read_cube,
    "read_cif": read_cif,
    "read_orca": read_orca,
    "read_charges": read_charges,
    "build_cluster": build_cluster,
    "print_cluster_xyz": reader.cif.Cluster.print_cluster_xyz,
    "print_cluster_kujo": reader.cif.Cluster.print_cluster_readable,
    "set_hard_cutoff": set_hard_cutoff,
    "set_int_cutoff": set_int_cutoff,
    "set_max_workers": set_max_workers,
    "generate_disorder": generate_disorder
}

if __name__ == "__main__":
    execute_time = time.time()
    instructions, options = utility.kujo_io.read_input("input.txt")
    suffix = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    out = f"/output{suffix}.out"
    for i in range(len(instructions)):
        str_execute_time = time.time()
        if len(options) != 0:
            result = dispatcher[instructions[i]](options[i])
        else:
            result = dispatcher[instructions[i]]()
        full_path = getcwd() + out
        finish_time = round((time.time() - str_execute_time) / 60, 4)
        if result is not None:
            file = open(full_path, "a+")
            file.write(f"{instructions[i]} with options {options[i]} returned ({finish_time} mins): {result}\n")
            file.close()
        else:
            file = open(full_path, "a+")
            file.write(f"Executed {instructions[i]} with options {options[i]}\n")
            file.close()
    overall_time = round((time.time() - execute_time) / 60, 4)
    file = open(full_path, "a+")
    file.write(f"Finished after {overall_time} minutes")
    file.close()
