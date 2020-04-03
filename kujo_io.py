from os import getcwd
import reader_cif


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
        words = x.split("=")
        if len(words) == 1:
            output_error("Syntax error at the input file: " + x, -1)
        instructions.append(words[0])
        options.append(words[1].split(","))
    return instructions, options


def print_cluster_xyz(cl: reader_cif.Cluster, path: str):
    full_path = getcwd()
    file = open(full_path + "/cluster.xyz", "w")
    for n in range(len(cl.molecules)):
        for m in range(cl.molecules[n].num_atoms):
            l = cl.molecules[n].atom_label[m]
            x = repr(round(cl.molecules[n].atom_coord[m, 0], 6))
            y = repr(round(cl.molecules[n].atom_coord[m, 1], 6))
            z = repr(round(cl.molecules[n].atom_coord[m, 2], 6))
            file.write(l + x + y + z)
    file.close()


