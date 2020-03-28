import excitons as e
import reader_cif as rc
import reader_orca as ro
import numpy as np
import linker

orca_molecule, mu = ro.read_orca("D:\[work]\PyGron\kes48_1.out")
cl = rc.Cluster(1, 1, 1, "D:\[work]\Kujo\KES48.cif")
mu_list = linker.link(cl, orca_molecule, mu)
#mu = np.zeros((1, 3))
#mu[0, 0] = 5.09777
#mu[0, 1] = 0.07312
#mu[0, 2] = 0.04931
#H = np.zeros((len(cl.mass_centers), len(cl.mass_centers)))
#d = 16.2516 # 8.6 angstrom
#e.hamiltonian_extended_dipole(cl, d, mu, H)
#e.spectra(cl, mu, H, 1000, 2000)
