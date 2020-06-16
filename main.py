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
                    i = int(opt[i1 + 1])
                    dsp[opt[i1]] = i
                except ValueError:
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


def read_cube(options_dispatcher: dict):
    global cube
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cube = reader.cube.Cube(full_path)


def read_cif(options_dispatcher: dict):
    global cif
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    cif = reader.cif.CifFile(full_path)


def read_charges(options_dispatcher: dict):
    global charges
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    charges = reader.charges.Charges(full_path)


def read_orca(options_dispatcher: dict):
    global orca
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    orca = reader.orca.Orca(full_path)


def cube_integration_wrapper(par_list: list):
    if max_w is not None:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_w) as executor:
            results = executor.map(reader.cube.integrate_cubes, par_list)
            for x in results:
                r = r + x
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(reader.cube.integrate_cubes, par_list)
            for x in results:
                r = r + x
    return r


def translated_coupling_td_integration(options_dispatcher: dict):
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
    mol1 = deepcopy(cube.molecule)
    mol2 = deepcopy(cube.molecule)
    mol2.atom_coord = mol2.atom_coord + translate
    for i1 in range(cube.steps[0, 0]):
        it.append([mol1, mol2, cube, i1])
    return cube_integration_wrapper(it)


def set_hard_cutoff(options_dispatcher: dict):
    global hard_cutoff
    hard_cutoff = float(options_dispatcher["value"])


def set_int_cutoff(options_dispatcher: dict):
    global int_cutoff
    int_cutoff = float(options_dispatcher["value"])


def set_max_workers(options_dispatcher: dict):
    global max_w
    max_w = int(options_dispatcher["value"])


def build_cluster(options_dispatcher: dict):
    global cluster
    global cif
    cluster = reader.cif.Cluster(cif, options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])
    cluster.build()
    it = []
    cluster.to_cartesian()
    if not ("print_cluster_xyz_premolecules" in instructions):
        cluster.connectivity()
        cluster.rebuild()
        cluster.find_rotations()
        cluster.multiply(options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])
        cluster.build_rmc()


def start_dipole_extended_dipole(options_dispatcher: dict):
    options_dispatcher_temp = deepcopy(options_dispatcher)
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
    if options_dispatcher["mu_x"] != 0.0 and options_dispatcher["mu_y"] != 0.0 and options_dispatcher["mu_z"] != 0.0:
        mu = np.array([options_dispatcher["mu_x"], options_dispatcher["mu_y"], options_dispatcher["mu_z"]])
    else:
        mu = orca.mu
    return t, mu


def translated_coupling_extended_dipole(options_dispatcher: dict):
    t, mu, d = start_dipole_extended_dipole(options_dispatcher)
    d = options_dispatcher["d"]
    mol1 = cluster.molecules[0]
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    return energy.excitons.coupling_extended_dipole(mol1, mol2, mu, d)


def translated_coupling_dipole(options_dispatcher: dict):
    t, mu = start_dipole_extended_dipole(options_dispatcher)
    mol1 = cluster.molecules[0]
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    return energy.excitons.coupling_dipole(mol1, mol2, mu)


def translated_coupling_charges(options_dispatcher: dict):
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


def generate_disorder(options_dispatcher: dict):
    global disorders
    n = H.shape[0]
    for i in range(options_dispatcher["n"]):
        disorders.append(energy.excitons.diagonal_disorder_sample(n, options_dispatcher["sigma"]))


def calculate_coupling(options_dispatcher: dict):
    if options_dispatcher["mol1"] is None and options_dispatcher["mol2"] is None:
        mol1 = cluster.molecules[options_dispatcher["site1"]]
        mol2 = cluster.molecules[options_dispatcher["site2"]]
    else:
        mol1 = options_dispatcher["mol1"]
        mol2 = options_dispatcher["mol2"]
    if options_dispatcher["mu_x"] == 0.0 and options_dispatcher["mu_y"] == 0.0 and options_dispatcher["mu_z"]:
        if options_dispatcher["method"] != "charges" and options_dispatcher["method"] != "integration":
            mu = orca.mu
    else:
        mu = np.array([options_dispatcher["mu_x"], options_dispatcher["mu_y"], options_dispatcher["mu_z"]])
    if options_dispatcher["method"] == "dipole":
        return energy.excitons.coupling_dipole(mol1, mol2, mu)
    elif options_dispatcher["method"] == "extended_dipole":
        return energy.excitons.coupling_extended_dipole(mol1, mol2, mu, options_dispatcher["d"])
    elif options_dispatcher["method"] == "charges":
        return energy.excitons.coupling_charges(mol1, mol2, charges.q)
    elif options_dispatcher["method"] == "integration":
        it = []
        for i1 in range(cube.steps[0, 0]):
            it.append([mol1, mol2, cube, i1])
        return cube_integration_wrapper(it)


def calculate_hamiltonian(options_dispatcher: dict):
    global H
    if options_dispatcher["periodic"]:
        cluster.build_rmc_periodic()
    else:
        cluster.build_rmc()
    for n in range(len(cluster.molecules)):
        for m in range(n + 1, len(cluster.molecules)):
            mol1 = deepcopy(cluster.molecules(n))
            mol2 = deepcopy(mol1)
            mol2.atom_coord += cluster.r_matrix[n][m]
            reader.cif.transform(mol2, mol2.rotation)
            v_out = v
            v_out.append("mol1")
            v_out.append(mol1)
            v_out.append("mol2")
            v_out.append(mol2)
            H[n, m] = calculate_coupling(v_out)
            H[m, n] = H[n, m]


def print_dimer_wrapper(options_dispatcher: dict):
    utility.kujo_io.print_dimer(cluster, options_dispatcher["site1"], options_dispatcher["site2"])


def print_cluster_xyz(options_dispatcher: dict):
    full_path = getcwd()
    file = open(full_path + "/cluster.xyz", "w")
    total_number = 0
    for n in range(len(cluster.molecules)):
        total_number += cluster.molecules[n].num_atoms
    file.write(repr(total_number) + "\n")
    file.write("XYZ file of molecular cluster generated in Kujo\n")
    for n in range(len(cluster.molecules)):
        for m in range(cluster.molecules[n].num_atoms):
            l = cluster.molecules[n].atom_label[m]
            x = repr(round(cluster.molecules[n].atom_coord[m, 0], 6))
            y = repr(round(cluster.molecules[n].atom_coord[m, 1], 6))
            z = repr(round(cluster.molecules[n].atom_coord[m, 2], 6))
            file.write(f'{l}   {x}   {y}   {z}\n')
    file.close()


def print_cluster_xyz_premolecules(options_dispatcher: dict):
    full_path = getcwd()
    file = open(full_path + "/cluster_premolecules.xyz", "w")
    total_number = 0
    for n in range(len(cluster.pre_molecules)):
        total_number += cluster.pre_molecules[n].num_atoms
    file.write(repr(total_number) + "\n")
    file.write("XYZ file of molecular cluster generated in Kujo\n")
    for n in range(len(cluster.pre_molecules)):
        for m in range(cluster.pre_molecules[n].num_atoms):
            l = cluster.pre_molecules[n].atom_label[m]
            x = repr(round(cluster.pre_molecules[n].atom_coord[m, 0], 6))
            y = repr(round(cluster.pre_molecules[n].atom_coord[m, 1], 6))
            z = repr(round(cluster.pre_molecules[n].atom_coord[m, 2], 6))
            file.write(f'{l}   {x}   {y}   {z}\n')
    file.close()


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
    "print_cluster_kujo": reader.cif.Cluster.print_cluster_readable,
    "set_hard_cutoff": set_hard_cutoff,
    "set_int_cutoff": set_int_cutoff,
    "set_max_workers": set_max_workers,
    "generate_disorder": generate_disorder,
    "calculate_coupling": calculate_coupling,
    "print_dimer": utility.kujo_io.print_dimer,
    "print_cluster_xyz": print_cluster_xyz,
    "print_cluster_xyz_premolecules": print_cluster_xyz_premolecules
}

options_list = {
    "read_cube": ["file"],
    "read_cif": ["file"],
    "read_charges": ["file"],
    "read_orca": ["file"],
    "translated_coupling_td_integration": ["vector", "vector_cif", "multiplier"],
    "set_hard_cutoff": ["value"],
    "set_int_cutoff": ["value"],
    "set_max_workers": ["value"],
    "build_cluster": ["a", "b", "c"],
    "translated_coupling_extended_dipole": ["vector", "vector_cif", "multiplier", "d", "mu_x", "mu_y", "mu_z"],
    "translated_coupling_dipole": ["vector", "vector_cif", "multiplier", "d", "mu_x", "mu_y", "mu_z"],
    "translated_coupling_charges": ["vector", "vector_cif", "multiplier"],
    "generate_disorder": ["sigma", "n"],
    "calculate_coupling": ["method", "site1", "site2", "mu_x", "mu_y", "mu_z", "d", "mo11", "mol2"],
    "calculate_hamiltonian": ["method", "periodic", "mu_x", "mu_y", "mu_z", "d"],
    "print_dimer_wrapper": ["site1", "site2"],
    "print_cluster_xyz": [],
    "print_cluster_xyz_premolecules": []
}


if __name__ == "__main__":
    execute_time = time.time()
    instructions, options = utility.kujo_io.read_input("input.txt")
    suffix = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    out = f"/output{suffix}.out"
    for i in range(len(instructions)):
        str_execute_time = time.time()
        if len(options) != 0:
            opt_to_method = {}
            l = options_list[instructions[i]]
            for i1 in range(len(l)):
                opt_to_method[l[i]] = ""
            options_parse(opt_to_method, options[i])
            result = dispatcher[instructions[i]](opt_to_method)
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
