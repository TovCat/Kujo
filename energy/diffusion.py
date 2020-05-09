import numpy as np
import scipy.constants


A = 219500  # Hartree to cm-1
bohr = 0.529177210  # bohr in Å


def participation_ratio(H: np.array, E: np.array, c: np.array):
    numerators = []
    denominators = []
    for q1 in range(H.shape[0]):
        for q2 in range(q1 + 1, H.shape[0]):
            inner = 0
            for n in range(c[q2].shape[0]):
                inner += c[q2, n] ** 4
            if q1 != q2:
                numerators.append((E[q1] - E[q2]) * inner)
                denominators.append(E[q1] - E[q2])
    numerators_sum = 0
    denominators_sum = 0
    for i in range(len(numerators)):
        numerators_sum += numerators[i]
        denominators_sum += denominators[i]
    return 1 / ((numerators_sum / len(numerators)) / (denominators_sum / len(denominators)))


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
    theta_linspace = np.linspace(0, np.pi, bins)
    theta_weights = np.zeros(bins)
    theta_sum = 0.0
    for i in bins:
        theta_weights[i] = np.sin(theta_linspace[i])
        theta_sum += theta_weights[i]
    theta_weights /= theta_sum
    theta = np.random.choice(theta_linspace, bins, p=theta_weights)
    phi = np.random.uniform(0, 2 * np.pi, bins)
    distribution = []
    for x in bins:
        for y in bins:
            theta_sin = np.sin(theta[x])
            theta_cos = np.cos(theta[x])
            phi_sin = np.sin(phi[y])
            phi_cos = np.cos(phi[y])
            distribution.append(np.array([theta_sin * phi_cos, theta_sin * phi_sin, theta_cos]))
    return distribution


def diffusion_no_thermal(H: np.array, E: np.array, c: np.array, r: np.array, u: np.array, gamma):
    D = 0.0
    for n in range(H.shape[0]):
        for m in range(n + 1, H.shape[0]):
            for n_j in range(H.shape[0]):
                for m_j in range(H.shape[0]):
                    j = np.inner(u, r[n, m] / bohr) * (H[n_j, m_j] / A) * c[n, n_j] * c[m, m_j]
                    D += (gamma / (gamma ** 2 + (E[n] / A - E[m] / A) ** 2)) * (j ** 2)
    return D / H.shape[0]


def diffusion_thermal(H: np.array, E: np.array, c: np.array, r: np.array, u: np.array, gamma, T):
    """
    Haken–Strobl–Reineker model
    """
    Z = 0.0
    D = 0.0
    exp_sum = 0.0
    for n in range(H.shape[0]):
        exponent = np.exp((-1 * scipy.constants.hbar * E[n]) / (scipy.constants.Boltzmann * T))
        exp_sum += exponent
        for m in range(n + 1, H.shape[0]):
            for n_j in range(H.shape[0]):
                for m_j in range(H.shape[0]):
                    pre_coeff = gamma / (gamma ** 2 + (E[n] - E[m]) ** 2)
                    j = np.inner(u, r[n, m] / bohr) * (H[n_j, m_j] / A) * c[n, n_j] * c[m, m_j]
                    D += pre_coeff * j * j * exponent
        return D / exp_sum
