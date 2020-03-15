import excitons as e
import reader_cif as rc
import reader_orca as ro
import numpy as np

orca_path = input("ORCA file path: ")
ro.read_orca(orca_path) # "D:\[work]\PyGron\kes48_1.out"
file_path = input("CIF file path: ")
cl = rc.Cluster(1, 1, 1, file_path)  # D:\[work]\CIF\p21c\FP4.cif D:\[work]\Kujo\KES48.cif
mu = np.zeros((1, 3))
mu[0, 0] = 5.09777
mu[0, 1] = 0.07312
mu[0, 2] = 0.04931
H = np.zeros((len(cl.molecules), len(cl.molecules)))
e.hamiltonian_dipole(cl, mu, H)
e.spectra(cl, mu, H)
