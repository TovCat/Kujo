from os import getcwd
import reader.cif
import reader.cube
import time


inp = ""
out = ""
version = "24.07.alpha"
parameters = {
}


def clear_string(a: str):
    return a.replace(" ", "").replace("\t", "").replace("\n", "")


def output_error(text: str, error_code: int):
    file = open(out, "w")
    crash_time = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    error_output = f"""
        !!!CALCULATION CRASHED!!!
        Reason: {text}
        Time: {crash_time}
        """
    file.write(error_output)
    file.close()
    exit(error_code)


def output_result(instruction: str, options: list, result, finish_time):
    file = open(out, "w")
    if result is not None:
        result_output = f"""
            Calculation finished!
            Routine: {instruction}
            Options: """
    else:
        result_output = f"""Routine finished!
            Routine: {instruction}
            Options: """
    for i in range(len(options) // 2):
        result_output += f"\n{options[2 * i]} = {options[2 * i + 1]}"
    exec_time = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    if result is not None:
        result_output += f"\nResult: {result}\nExecution (wall) time: {finish_time} min.\nCurrent time and date: {exec_time}"
    else:
        result_output += f"\nExecution (wall) time: {finish_time} min.\nCurrent time and date: {exec_time}"
    file.write(result_output)
    file.close()


def initiate_output():
    file = open(out, "w")
    exec_time = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    for_write = f"""
        ====================================================================
        Kujo: Python code for exciton dynamics in organic single crystals
        
        Author: Igor Koskin
        Affiliation: Novosibirsk Institute of Organic Chemistry
        e-mail: osingran@yandex.ru
        
        Version: {version}
        ====================================================================

        Executed input file: {inp}
        Time and date: {exec_time}
        !!!DO NOT DELETE THE INPUT FILE. KUJO DOES NOT SAVE INPUT INSTRUCTIONS IN THE OUTPUT!!!

        """
    file.write(for_write)
    file.close()


def finalize_output(overall_time):
    file = open(out, "w")
    exec_time = f"-{str(time.localtime().tm_mday)}-{str(time.localtime().tm_mon)}-{str(time.localtime().tm_year)}-{str(time.localtime().tm_hour)}-{str(time.localtime().tm_min)}-{str(time.localtime().tm_sec)}"
    for_write = f"""
        ====================================================================
        CALCULATION TERMINATED NORMALLY

        Input file: {inp}
        Overall execution wall time: {overall_time}
        Time and date: {exec_time}
        ====================================================================
        """
    file.write(for_write)
    file.close()


def read_options():
    global parameters
    contents = []
    try:
        full_path = getcwd() + "/config.ini"
        file = open(full_path, "r")
        contents = file.readlines()
        file.close()
    except IOError:
        exit(-1)
    for i in range(len(contents)):
        words = contents[i].split("=")
        try:
            parameters[words[0]] = float(words[1])
        except ValueError:
            exit(-1)


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
        output_error("Empty input file.", -1)
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