from os import getcwd
import reader.cif
import reader.cube


def clear_string(a: str):
    return a.replace(" ", "").replace("\t", "").replace("\n", "")


def output_error(text: str, error_code: int):
    full_path = getcwd()
    file = open(full_path + "/error_msg.txt", "w")
    file.write(text)
    file.close()
    exit(error_code)


def read_input(path: str):
    contents = []
    try:
        file = open(path, "r")
        contents = file.readlines()
        file.close()
    except OSError:
        output_error("Could not open input file at: " + path, -1)
    instructions = []
    options = []
    if len(contents) == 0:
        output_error("Empty input file!", -1)
    for_pop = []
    for x in contents:
        contents[contents.index(x)] = clear_string(x).lower()
        if contents[contents.index(x)] == "" or contents[contents.index(x)][0] == "#":
            for_pop.append(contents.index(x))
    for i in range(len(for_pop)):
        contents.pop(for_pop[i])
    for_slice = [0]
    sliced = []
    for i in range(len(contents)):
        if contents[i] == "end":
            for_slice.append(i)
            temp = []
            for i1 in range(for_slice[len(for_slice) - 2], for_slice[len(for_slice) - 1]):
                temp.append(contents[i])
            sliced.append(temp)
    for i in range(len(sliced)):
        instructions.append(sliced[i][0])
        for i1 in range(len(sliced[i])):
            words = sliced[i][i1].split("=")
            try:
                a = float(words[1])
            except ValueError:
                if "," in words[1]:
                    a = words[1].split(",")
                    for x in range(len(a)):
                        try:
                            a[x] = float(a[x])
                        except ValueError:
                            pass
            options.append([words[0], a])
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