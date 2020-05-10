import numpy as np
import reader.cif


class Orca:

    def __init__(self, path=""):
        file = open(path, "r")
        contents = file.readlines()
        file.close()
        index = -1
        self.mu = np.zeros((1, 3))
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
                index = contents.index(x)
                s = contents[index + 5]
                words = s.split()
                for i in range(3):
                    self.mu[0, i] = float(words[i + 5])
                break
        self.molecule = reader.cif.Molecule(len(coord_list))
        for x in coord_list:
            index = coord_list.index(x)
            words = x.split()
            self.molecule.atom_label.append(words[0])
            for i in range(3):
                self.molecule.atom_coord[index, i] = float(words[i + 1])