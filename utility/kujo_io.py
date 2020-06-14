from os import getcwd
import reader.cif
import reader.cube


def output_error(text: str, error_code: int):
    full_path = getcwd()
    file = open(full_path + "/error_msg.txt", "w")
    file.write(text)
    file.close()
    exit(error_code)


def read_input(path: str):
    full_path = getcwd() + "/" + path
    contents = []
    try:
        file = open(full_path, "r")
        contents = file.readlines()
        file.close()
    except OSError:
        output_error("Could not open input file at: " + full_path, -1)
    instructions = []
    options = []
    if len(contents) == 0:
        output_error("Empty input file!", -1)
    for x in contents:
        if ":" in x:
            words = x.split(":")
            instructions.append(words[0].strip())
            words[1] = words[1].strip()
            pre_options1 = words[1].split(";")
            flag = 0
            for x in pre_options1:
                x = x.strip()
                if len(pre_options1) == 1:
                    options.append(x.split("="))
                elif flag == 0:
                    options.append(x.split("="))
                    flag = len(options) - 1
                else:
                    pre_options2 = x.split("=")
                    for y in pre_options2:
                        options[flag].append(y)
        else:
            instructions.append(x)
    if len(options) != 0:
        for x in range(len(options)):
            for y in range(len(options[x])):
                options[x][y] = options[x][y].strip()
    return instructions, options


def print_hamiltonian(H):
    full_path = getcwd()
    file = open(full_path + "/hamiltonian.kujo", "w")
    for n in range(H.shape[0]):
        for m in range(H.shape[1]):
            file.write(repr(round(H[n, m], 6)) + " ")
        file.write("\n")
    file.close()


def print_dimer(cl: reader.cif.Cluster, s1, s2: int):
    full_path = getcwd()
    file = open(f'{full_path}/site_{s1}_{s2}.xyz', "w")
    total_number = 0
    for n in range(len(cl.molecules)):
        total_number += + cl.molecules[n].num_atoms
    file.write(repr(total_number) + "\n")
    file.write("XYZ file of molecular cluster generated in Kujo\n")
    for n in range(len(cl.molecules)):
        for m in range(cl.molecules[n].num_atoms):
            if n != s1 or n != s2:
                l = "H"
            else:
                l = cl.molecules[n].atom_label
            x = repr(round(cl.molecules[n].atom_coord[m, 0], 6))
            y = repr(round(cl.molecules[n].atom_coord[m, 1], 6))
            z = repr(round(cl.molecules[n].atom_coord[m, 2], 6))
            file.write(f'{l}   {x}   {y}   {z}\n')
    file.close()