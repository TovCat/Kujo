import numpy as np
import reader.cif
import reader.charges
from copy import deepcopy

A = 219500  # Hartree to cm-1
bohr = 0.529177210  # bohr in Å


def align_dipole(angles: np.array, mu: np.array):
    x_r, y_r, z_r = reader.cif.rotation_matrix(angles[0], angles[1], angles[2])
    temp = np.transpose(np.matmul(np.transpose(mu), x_r))
    temp = np.transpose(np.matmul(np.transpose(temp), y_r))
    temp = np.transpose(np.matmul(np.transpose(temp), z_r))
    return temp


def coupling_dipole(mol1: reader.cif.Molecule, mol2: reader.cif.Molecule, mu: np.array,):
    mu_1 = mu
    mu_2 = mu
    reader.cif.transform(mu_1, mol1.rotation)
    reader.cif.transform(mu_2, mol2.rotation)
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
    reader.cif.transform(mu_1, mol1.rotation)
    reader.cif.transform(mu_2, mol2.rotation)
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


def coupling_charges(c: reader.charges.Charges, mol_temp: reader.cif.Molecule, dr, r):
    mol = deepcopy(c.mol)
    x_rot, y_rot, z_rot = reader.cif.rotation_matrix(mol_temp.alpha, mol_temp.beta, mol_temp.gamma)
    reader.cif.transform(mol, x_rot)
    reader.cif.transform(mol, y_rot)
    reader.cif.transform(mol, z_rot)
    mol.atom_coord += dr
    J = 0.0
    for n in range(c.mol.num_atoms):
        for m in range(c.mol.num_atoms):
            J += c.q[n] * c.q[m] / r
    return J * A


def hamiltonian_diagonal_disorder(H, sigma, exp):
    for n in range(H.shape[0]):
        H[n, n] = np.random.normal(exp, sigma)