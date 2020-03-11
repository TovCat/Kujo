import numpy as np
import reader_cif as rc


def read_orca(path: str):
    file = open(path, "r")
    contents = file.readlines()
    file.close()
    index = -1
    mu = np.zeros((3, 1))
    for x in contents:
        index = index + 1
        x.strip()
        if x.find("CARTESIAN COORDINATES (ANGSTROEM)") != -1:
            coord_list = []
            index = index + 2
            s = contents[index]
            while s != '\n':
                s.strip()
                coord_list.append(s)
                index = index + 1
                s = contents[index]
        if x.find("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS") != -1:
            s = contents[index + 5]
            words = s.split()
            for i in range(3):
                mu[i, 0] = float(words[i + 5])
                break
    orca_molecule = rc.Molecule(len(coord_list))
    for x in coord_list:
        index = coord_list.index(x)
        words = x.split()
        orca_molecule.atom_label.append(words[0])
        for i in range(3):
            orca_molecule.atom_coord[index, i] = float(words[i + 1])

    return orca_molecule, mu
