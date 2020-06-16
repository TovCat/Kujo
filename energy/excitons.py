import numpy as np
import reader.cif
import reader.charges
from copy import deepcopy

A = 219500  # Hartree to cm-1
bohr = 0.529177210  # bohr in Å


def coupling_dipole(mol1: reader.cif.Molecule, mol2: reader.cif.Molecule, mu: np.array):
    mu_1 = mu
    mu_2 = mu
    mu_1 = np.transpose(np.matmul(np.transpose(mu_1), mol1.rotation))
    mu_2 = np.transpose(np.matmul(np.transpose(mu_2), mol2.rotation))
    r = mol1.mass_center() - mol2.mass_center()
    r = r / bohr
    r1 = np.linalg.norm(r) / bohr
    r3 = r1 ** 3
    r5 = r1 ** 5
    J = np.inner(mu_1, mu_2) / r3 - 3 * (np.inner(mu, r) * np.inner(r, mu)) / r5
    return J * A


def coupling_extended_dipole(mol1: reader.cif.Molecule, mol2: reader.cif.Molecule, mu: np.array, d):
    mu_1 = mu
    mu_2 = mu
    mu_1 = np.transpose(np.matmul(np.transpose(mu_1), mol1.rotation))
    mu_2 = np.transpose(np.matmul(np.transpose(mu_2), mol2.rotation))
    d = d / bohr
    q = np.linalg.norm(mu_1) / (2 * d)
    mu_trans_1 = mu_1 * (d / np.linalg.norm(mu_1))
    mu_trans_2 = mu_2 * (d / np.linalg.norm(mu_2))
    r = mol1.mass_center() - mol2.mass_center()
    r_bohr = r / bohr
    p1_1 = mu_trans_1
    p1_2 = -1 * mu_trans_1
    p2_1 = r_bohr + mu_trans_2
    p2_2 = r_bohr - mu_trans_2
    r_pp = np.linalg.norm(p1_1 - p2_1)
    r_pm = np.linalg.norm(p1_1 - p2_2)
    r_mp = np.linalg.norm(p1_2 - p2_1)
    r_mm = np.linalg.norm(p1_2 - p2_2)
    return (q ** 2) * ((1 / r_pp) - (1 / r_pm) - (1 / r_mp) + (1 / r_mm)) * A


def coupling_charges(mol1: reader.cif.Molecule, mol2: reader.cif.Molecule, q):
    reader.cif.transform(mol1.atom_coord, mol1.rotation)
    reader.cif.transform(mol2.atom_coord, mol2.rotation)
    J = 0.0
    for n in range(mol1.num_atoms):
        for m in range(mol2.num_atoms):
            r = np.linalg.norm(mol1.atom_coord[n:n + 1] - mol2.atom_coord[m:m + 1]) / bohr
            J += q[n] * q[m] / r
    return J * A


def diagonal_disorder_sample(n, sigma):
    l = []
    for i in range(n):
        l.append(np.random.normal(0, sigma))
    return l