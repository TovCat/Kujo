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
import energy.diffusion

cube = reader.cube.Cube()
cif = reader.cif.CifFile()
cluster = reader.cif.Cluster()
charges = reader.charges.Charges()
orca = reader.orca.Orca()
H = np.zeros((3, 3))
disorders = []
parameters = {
}


def read_options():
    global parameters
    full_path = ""
    contents = []
    try:
        full_path = getcwd() + "/config.ini"
        file = open(full_path, "r")
        contents = file.readlines()
        file.close()
    except IOError:
        utility.kujo_io.output_error("Could not open input file at: " + full_path, -1)
    for i in range(len(contents)):
        words = contents[i].split("=")
        try:
            parameters[words[0]] = float(words[1])
        except ValueError:
            exit(-1)


def options_parse(dsp: dict, opt: list):
    for i_op in range(0, len(opt), 2):
        if opt[i_op] in dsp:
            if opt[i_op + 1].lower() == "true":
                dsp[opt[i_op]] = True
            elif opt[i_op + 1].lower() == "false":
                dsp[opt[i_op]] = False
            else:
                try:
                    i_test = int(opt[i_op + 1])
                    dsp[opt[i_op]] = i_test
                except ValueError:
                    try:
                        fl = float(opt[i_op + 1])
                        dsp[opt[i_op]] = fl
                    except ValueError:
                        if "[" in opt[i_op + 1] or "]" in opt[i_op + 1]:
                            opt[i_op + 1].replace("[", "")
                            opt[i_op + 1].replace("]", "")
                            words = opt[i_op + 1].split(",")
                            temp = np.zeros([len(words), 1])
                            try:
                                for i_op2 in range(len(words)):
                                    temp[i_op2, 0] = float(words[i])
                                dsp[opt[i_op]] = temp
                            except ValueError:
                                exit(-1)
                        else:
                            dsp[opt[i_op]] = opt[i_op + 1]
        else:
            exit(-1)


def read_file(options_dispatcher: dict):
    global cube
    global cif
    global charges
    global orca
    full_path = ""
    file_types = {
        "cube": cube.read,
        "cub": cube.read,
        "cif": cif.read,
        "chg": charges.read,
        "out": orca.read
    }
    w = options_dispatcher['file'].split('.')
    full_path = getcwd() + f"/{options_dispatcher['file']}"
    file_types[w[1]](full_path)


def cube_integration_wrapper(par_list: list):
    if parameters["max_w"] != -1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["max_w"]) as executor:
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
        translate = 0
        exit(-1)
    it = []
    mol1 = deepcopy(cube.molecule)
    mol2 = deepcopy(cube.molecule)
    mol2.atom_coord = mol2.atom_coord + translate
    for i_c in range(cube.steps[0, 0]):
        it.append([mol1, mol2, cube, i_c])
    return cube_integration_wrapper(it)


def build_cluster(options_dispatcher: dict):
    global cluster
    global cif
    cluster.init_cif(cif)
    cluster.build()
    cluster.to_cartesian()
    if not ("print_cluster_xyz_premolecules" in instructions):
        cluster.connectivity()
        cluster.rebuild()
        cluster.find_rotations()
        cluster.multiply(options_dispatcher["a"], options_dispatcher["b"], options_dispatcher["c"])
        cluster.build_rmc()


def start_dipole_extended_dipole(options_dispatcher: dict):
    if options_dispatcher["vector_cif"] == "a":
        t = cif.vector_a * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "b":
        t = cif.vector_b * options_dispatcher["multiplier"]
    elif options_dispatcher["vector_cif"] == "c":
        t = cif.vector_c * options_dispatcher["multiplier"]
    elif options_dispatcher["vector"] is not None:
        t = options_dispatcher["vector"] * options_dispatcher["multiplier"]
    else:
        t = 0
        exit(-1)
    if options_dispatcher["mu_x"] != 0.0 and options_dispatcher["mu_y"] != 0.0 and options_dispatcher["mu_z"] != 0.0:
        mu = np.array([options_dispatcher["mu_x"], options_dispatcher["mu_y"], options_dispatcher["mu_z"]])
    else:
        mu = orca.mu
    return t, mu


def translated_coupling_extended_dipole(options_dispatcher: dict):
    t, mu = start_dipole_extended_dipole(options_dispatcher)
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
        t = 0
        exit(-1)
    mol1 = charges.mol
    mol2 = deepcopy(mol1)
    mol2.atom_coord = mol2.atom_coord + t
    return energy.excitons.coupling_charges(mol1, mol2, charges.q)


def generate_disorder(options_dispatcher: dict):
    global disorders
    disorders = []
    n = len(cluster.molecules)
    for i_d in range(options_dispatcher["n"]):
        disorders.append(energy.excitons.diagonal_disorder_sample(n, options_dispatcher["sigma"]))


def calculate_coupling(options_dispatcher: dict):
    if options_dispatcher["mol1"] is None and options_dispatcher["mol2"] is None:
        mol1 = cluster.molecules[options_dispatcher["site1"]]
        mol2 = cluster.molecules[options_dispatcher["site2"]]
    else:
        mol1 = options_dispatcher["mol1"]
        mol2 = options_dispatcher["mol2"]
    mu = np.zeros((3))
    if options_dispatcher["mu_x"] == 0.0 and options_dispatcher["mu_y"] == 0.0 and options_dispatcher["mu_z"] \
            and options_dispatcher["method"] != "charges" and options_dispatcher["method"] != "integration":
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
    H = np.zeros((len(cluster.molecules), len(cluster.molecules)))
    for n in range(len(cluster.molecules)):
        for m in range(n + 1, len(cluster.molecules)):
            mol1 = deepcopy(cluster.molecules[n])
            mol2 = deepcopy(cluster.molecules[m])
            dict_out = options_dispatcher
            dict_out["mol1"] = mol1
            dict_out["mol2"] = mol2
            H[n, m] = calculate_coupling(dict_out)
            H[m, n] = H[n, m]


def print_dimer_wrapper(options_dispatcher: dict):
    utility.kujo_io.print_dimer(cluster, options_dispatcher["site1"], options_dispatcher["site2"])


def print_cluster_xyz(options_dispatcher: dict):
    def print_cluster(mol: reader.cif.Molecule, f):
        for m in range(mol.num_atoms):
            l = mol.atom_label[m]
            x = repr(round(mol.atom_coord[m, 0], 6))
            y = repr(round(mol.atom_coord[m, 1], 6))
            z = repr(round(mol.atom_coord[m, 2], 6))
            f.write(f'{l}   {x}   {y}   {z}\n')
    
    full_path_xyz = getcwd()
    file_xyz = open(full_path_xyz + "/cluster.xyz", "w")
    if options_dispatcher["mode"] == "premolecules":
        total_number_xyz = 0
        for n in range(len(cluster.pre_molecules)):
            total_number_xyz += cluster.pre_molecules[n].num_atoms
    else:
        total_number_xyz = 0
        for n in range(len(cluster.molecules)):
            total_number_xyz += cluster.molecules[n].num_atoms
    file_xyz.write(repr(total_number_xyz) + "\n")
    file_xyz.write("XYZ file of molecular cluster generated in Kujo\n")
    if options_dispatcher["mode"] == "premolecules":
        for n in range(len(cluster.pre_molecules)):
            print_cluster(cluster.pre_molecules[n], file_xyz)
    else:
        for n in range(len(cluster.molecules)):
            print_cluster(cluster.molecules[n], file_xyz)
    file_xyz.close()


def calculate_participation_ratio(options_dispatcher: dict):
    p = []
    if len(disorders) != 0:
        for i in range(len(disorders)):
            for n in range(len(cluster.molecules)):
                H[n, n] = disorders[i][n]
            E, c = np.linalg.eig(H)
            p.append(energy.diffusion.participation_ratio(H, E, c))
    p_sum = 0
    for i_pr in range(len(p)):
        p_sum += p[i_pr]
    return p_sum / len(p)


def calculate_diffusion(options_dispatcher: dict):
    if options_dispatcher["mode"] == "xy" or options_dispatcher["mode"] == "yz" or options_dispatcher["mode"] == "xz":
        distribution = energy.diffusion.distribute_over_plane(options_dispatcher["mode"], options_dispatcher["bins"])
    elif options_dispatcher["mode"] == "sphere":
        distribution = energy.diffusion.distribute_over_sphere(options_dispatcher["bins"])
    elif options_dispatcher["mode"] == "a":
        distribution = [np.array([1, 0, 0])]
    elif options_dispatcher["mode"] == "b":
        distribution = [np.array([0, 1, 0])]
    elif options_dispatcher["mode"] == "a":
        distribution = [np.array([0, 0, 1])]
    else:
        distribution = []
        exit(-1)
    T = options_dispatcher["temperature"]
    gamma = options_dispatcher["gamma"]
    par_list = []
    diffusion_different_samples = []
    for i_outer in range(len(disorders)):
        H_disorder = H
        diffusion_specific_disorder = []
        for i_d in range(len(cluster.molecules)):
            H_disorder[i_d, i_d] = disorders[i_outer][i_d]
        E, c = np.linalg.eig(H_disorder)
        for i_cd in range(len(distribution)):
            if options_dispatcher["thermal"]:
                par_list.append([H_disorder, E, c, cluster.r_matrix, distribution[i_cd], gamma, T])
            else:
                par_list.append([H_disorder, E, c, cluster.r_matrix, distribution[i_cd], gamma])
        if parameters["max_w"] != -1:
            if options_dispatcher["thermal"]:
                with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["max_w"]) as executor:
                    results = executor.map(energy.diffusion.diffusion_thermal, par_list)
                    for x in results:
                        diffusion_specific_disorder.append(x)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["max_w"]) as executor:
                    results = executor.map(energy.diffusion.diffusion_no_thermal, par_list)
                    for x in results:
                        diffusion_specific_disorder.append(x)
        else:
            if options_dispatcher["thermal"]:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = executor.map(energy.diffusion.diffusion_thermal, par_list)
                    for x in results:
                        diffusion_specific_disorder.append(x)
            else:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = executor.map(energy.diffusion.diffusion_no_thermal, par_list)
                    for x in results:
                        diffusion_specific_disorder.append(x)
        diffusion_different_samples.append(diffusion_specific_disorder)
    diffusion_final = []
    for i_inner in range(options_dispatcher["bins"]):
        sum = 0
        for i_outer in range(len(disorders)):
            sum += diffusion_different_samples[i_outer][i_inner]
        diffusion_final.append(sum / len(disorders))
    return diffusion_final


dispatcher = {
    "translated_coupling_td_integration": translated_coupling_td_integration,
    "translated_coupling_extended_dipole": translated_coupling_extended_dipole,
    "translated_coupling_dipole": translated_coupling_dipole,
    "translated_coupling_charges": translated_coupling_charges,
    "read_file": read_file,
    "build_cluster": build_cluster,
    "print_cluster_kujo": reader.cif.Cluster.print_cluster_readable,
    "generate_disorder": generate_disorder,
    "calculate_coupling": calculate_coupling,
    "print_dimer": utility.kujo_io.print_dimer,
    "print_cluster_xyz": print_cluster_xyz,
    "calculate_participation_ratio": calculate_participation_ratio,
    "calculate_hamiltonian": calculate_hamiltonian,
    "calculate_diffusion": calculate_diffusion
}

options_list = {
    "read_file": ["file"],
    "translated_coupling_td_integration": ["vector", "vector_cif", "multiplier"],
    "build_cluster": ["a", "b", "c"],
    "translated_coupling_extended_dipole": ["vector", "vector_cif", "multiplier", "d", "mu_x", "mu_y", "mu_z"],
    "translated_coupling_dipole": ["vector", "vector_cif", "multiplier", "d", "mu_x", "mu_y", "mu_z"],
    "translated_coupling_charges": ["vector", "vector_cif", "multiplier"],
    "generate_disorder": ["sigma", "n"],
    "calculate_coupling": ["method", "site1", "site2", "mu_x", "mu_y", "mu_z", "d", "mo11", "mol2"],
    "calculate_hamiltonian": ["method", "periodic", "mu_x", "mu_y", "mu_z", "d"],
    "print_dimer_wrapper": ["site1", "site2"],
    "print_cluster_xyz": ["mode"],
    "calculate_participation_ratio": [],
    "calculate_diffusion": ["mode", "thermal", "bins", "gamma", "temperature"]
}


if __name__ == "__main__":
    execute_time = time.time()
    read_options()
    instructions, options = utility.kujo_io.read_input("input.txt")
    suffix = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    out = f"/output{suffix}.out"
    for i in range(len(instructions)):
        str_execute_time = time.time()
        if len(options) != 0:
            opt_to_method = {}
            l = options_list[instructions[i]]
            for i1 in range(len(l)):
                opt_to_method[l[i1]] = ""
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
