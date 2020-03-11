import excitons as e
import reader_cif as rc

file_path = input("CIF file path: ")
cl = rc.Cluster(1, 1, 1, file_path)  # D:\[work]\CIF\p21c\FP4.cif D:\[work]\Kujo\KES48.cif
e.hamiltonian(cl)
e.spectra(cl)
