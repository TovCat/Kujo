import numpy as np
import reader.cif

A = 219500  # Hartree to cm-1
bohr = 0.529177210  # bohr in Å


def align_dipole(angles: np.array, mu: np.array):
    x_r, y_r, z_r = reader.cif.rotation_matrix(angles[0], angles[1], angles[2])
    temp = np.transpose(np.matmul(np.transpose(mu), x_r))
    temp = np.transpose(np.matmul(np.transpose(temp), y_r))
    temp = np.transpose(np.matmul(np.transpose(temp), z_r))
    return temp


def coupling_dipole(angles1: np.array, angles2: np.array, mu: np.array, distance, r):
    mu_1 = align_dipole(angles1, mu)
    mu_2 = align_dipole(angles2, mu)
    r1 = distance / bohr
    r3 = r1 ** 3
    r5 = r1 ** 5
    J = np.inner(mu_1, mu_2) / r3 + (np.inner(mu, r) * np.inner(r, mu)) / r5
    return J * A


def coupling_extended_dipole(angles1: np.array, angles2: np.array, mu: np.array, r, d):
    mu_1 = align_dipole(angles1, mu)
    mu_2 = align_dipole(angles2, mu)
    d = d / bohr
    q = np.linalg.norm(mu) / (2 * d)
    mu_trans_1 = mu_1 * (d / np.linalg.norm(mu_1))
    mu_trans_2 = mu_2 * (d / np.linalg.norm(mu_2))
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
                u_v = np.array([0, np.sin(u[x]), np.cos(u[x])])
                r = (clust.mass_centers[n] - clust.mass_centers[m]) / bohr3
                j = np.inner(u_v, r) * H[n, m] / A
                D[x] = D[x] + (gamma / (gamma**2 + (E[n] / A - E[m] / A)**2)) * (j**2)
        D[x] = (1/len(clust.mass_centers)) * D[x]
        print(D[x])
