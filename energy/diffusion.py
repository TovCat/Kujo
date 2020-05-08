import numpy as np


A = 219500  # Hartree to cm-1
bohr = 0.529177210  # bohr in Å


def distribute_over_plane(plane, bins):
    u = np.linspace(0, 2 * np.pi, bins)
    distribution = []
    if plane == "xy":
        for x in range(bins):
            distribution.append(np.array([np.sin(u[x]), np.cos(u[x]), 0]))
    elif plane == "yz":
        for x in range(bins):
            distribution.append(np.array([0, np.sin(u[x]), np.cos(u[x])]))
    elif plane == "xz":
        for x in range(bins):
            distribution.append(np.array([np.sin(u[x]), 0, np.cos(u[x])]))
    return distribution


def distribute_over_sphere(bins):
    theta = np.random.uniform(0, np.pi, bins)
    phi = np.random.uniform(0, 2 * np.pi, bins)
    distribution = []
    for x in bins:
        theta_sin = np.sin(theta[x])
        theta_cos = np.cos(theta[x])
        phi_sin = np.sin(phi[x])
        phi_cos = np.cos(phi[x])
        distribution.append(np.array([theta_sin * phi_cos, theta_sin * phi_sin, theta_cos]))
    return distribution


def diffusion_no_thermal(H: np.array, E: np.array, r: np.array, u: np.array, gamma):
    D = 0.0
    for n in range(H.shape[0]):
        for m in range(n + 1, H.shape[0]):
            j = np.inner(u, r / bohr) * H[n, m] / A
            D += (gamma / (gamma ** 2 + (E[n] / A - E[m] / A) ** 2)) * (j ** 2)
    return D
