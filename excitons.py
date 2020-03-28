import matplotlib.pyplot as plt
import numpy as np
import time
import reader_cif as rc
import math

A = 219500  # Hartree to cm-1
bohr3 = 0.529177210  # bohr in Å

def hamiltonian_dipole(c: rc.Cluster, mu, H):
    print("Calculate Hamiltonian")
    h_time = time.time()
    for n in range(len(c.molecules)):
        for m in range(n + 1, len(c.molecules)):
            v1 = np.zeros((1, 3))
            v2 = np.zeros((1, 3))
            v1[0, 0] = c.mass_centers[n, 0]
            v1[0, 1] = c.mass_centers[n, 1]
            v1[0, 2] = c.mass_centers[n, 2]
            v2[0, 0] = c.mass_centers[m, 0]
            v2[0, 1] = c.mass_centers[m, 1]
            v2[0, 2] = c.mass_centers[m, 2]
            dv = v1 - v2
            r1 = np.linalg.norm(dv)
            r1 = r1 / bohr3
            r3 = r1 * r1 * r1
            r5 = r3 * r1 * r1
            J = np.inner(mu, mu) / r3 + (np.inner(mu, dv) * np.inner(dv, mu)) / r5
            H[n, m] = J * A
            H[m, n] = H[n, m]
    print("   Done: %s" % (time.time() - h_time))


def hamiltonian_extended_dipole(c: rc.Cluster, d, mu, H):
    q = np.linalg.norm(mu) / (2 * d)
    mu_trans = mu * math.sqrt(d / np.linalg.norm(mu))
    for n in range(len(c.mass_centers)):
        for m in range(n + 1, len(c.mass_centers)):
            mass_center_bohr1 = c.mass_centers[n:n + 1] / bohr3
            mass_center_bohr2 = c.mass_centers[m:m + 1] / bohr3
            p1_1 = mass_center_bohr1 + mu_trans
            p1_2 = mass_center_bohr1 - mu_trans
            p2_1 = mass_center_bohr2 + mu_trans
            p2_2 = mass_center_bohr2 - mu_trans
            r_pp = np.linalg.norm(p1_1 - p2_1)
            r_pm = np.linalg.norm(p1_1 - p2_2)
            r_mp = np.linalg.norm(p1_2 - p2_1)
            r_mm = np.linalg.norm(p1_2 - p2_2)
            J = (q ** 2) * ((1 / r_pp) - (1 / r_pm) - (1 / r_mp) + (1 / r_mm))
            H[n, m] = J * A
            H[m, n] = H[n, m]


def hamiltonian_diagonal_disorder(H, sigma, exp):
    for n in range(H.shape[0]):
        H[n, n] = np.random.normal(exp, sigma)


def spectra(clust: rc.Cluster, mu, H, bins, sigma):
    plt.xlabel('Energy (cm-1)')
    plt.ylabel('Intensity')
    mu_a = np.zeros((len(clust.mass_centers), 3))
    for x in range(len(clust.mass_centers)):
        mu_a[x, 0] = mu[0, 0]
        mu_a[x, 1] = mu[0, 1]
        mu_a[x, 2] = mu[0, 2]
    E, c = np.linalg.eig(H)
    EminA = np.min(E)
    EmaxA = np.max(E)
    print(EmaxA)
    Emin = EminA - 0.1 * (EmaxA - EminA) - 3 * sigma
    Emax = EmaxA + 0.1 * (EmaxA - EminA) + 3 * sigma
    dE = (Emax - Emin) / bins
    Ex = np.linspace(Emin, Emax, bins)
    Ey = np.zeros(bins)
    for n1 in range(len(clust.mass_centers)):
        bin = int(round((E[n1] - Emin) / dE))
        Emu = np.zeros(3)
        Emu[0] = np.inner(c[:, n1], mu_a[:, 0])
        Emu[1] = np.inner(c[:, n1], mu_a[:, 1])
        Emu[2] = np.inner(c[:, n1], mu_a[:, 2])
        Ey[bin] = Ey[bin] + np.linalg.norm(Emu) ** 2

    # Convolute spectrum
    # First create normalized Gaussian centered in the middle
    Cx = Ex - (Emax + Emin) / 2  # Make new axis with value zero in the middle of array
    Cy = np.exp(-Cx ** 2 / 2 / sigma ** 2) / np.sqrt(2 * np.pi * sigma ** 2)
    # Do the actual convolusion
    Ny = np.convolve(Ey, Cy, mode='same')
    # Plot everything in one final plot
    #plt.plot(Ex, Ey / np.max(Ey))
    # plt.plot(Ex,Cy/np.max(Cy))
    #plt.plot(Ex, Ny / np.max(Ny))
    #plt.show()

    u = np.linspace(0, 2*math.pi, bins)
    D = np.zeros(bins)
    gamma = (sigma / 2) / A
    for x in range(bins):
        for n in range(len(clust.mass_centers)):
            for m in range(n+1, len(clust.mass_centers)):
                u_v = np.array([0, math.sin(u[x]), math.cos(u[x])])
                r = (clust.mass_centers[n:n+1] - clust.mass_centers[m:m+1]) / bohr3
                j = np.inner(u_v, r) * H[n, m] / A
                D[x] = D[x] + (gamma / (gamma**2 + (E[n] / A - E[m] / A)**2)) * (j**2)
        D[x] = (1/len(clust.mass_centers)) * D[x]
        print(D[x])
