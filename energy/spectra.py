import numpy as np
import matplotlib.pyplot as plt


def common(E, c, mu, N, sigma):
    bins = 1000
    EminA = np.min(E)
    EmaxA = np.max(E)
    if (EmaxA-EminA) / sigma > 1000:
        bins = int(np.floor((EmaxA - EminA) / sigma))
    Emin = EminA - 0.1 * (EmaxA - EminA) - 3 * sigma
    Emax = EmaxA + 0.1 * (EmaxA - EminA) + 3 * sigma
    dE = (Emax - Emin) / bins
    Ex = np.linspace(Emin, Emax, bins)
    Ey = np.zeros(bins)
    for n in range(N):
        bin = int(round((E[n] - Emin) / dE))
        Emu = np.zeros(3)
        Emu[0] = np.inner(c[:, n], mu[:, 0])
        Emu[1] = np.inner(c[:, n], mu[:, 1])
        Emu[2] = np.inner(c[:, n], mu[:, 2])
        Ey[bin] = Ey[bin] + np.linalg.norm(Emu) ** 2
    Cx = Ex - (Emax + Emin) / 2
    Cy = np.exp(-Cx ** 2 / 2 / sigma ** 2) / np.sqrt(2 * np.pi * sigma ** 2)
    Ny = np.convolve(Ey, Cy, mode='same')
    plt.plot(Ex, Ey/np.max(Ey))
    plt.plot(Ex, Cy/np.max(Cy))
    plt.plot(Ex, Ny/np.max(Ny))
    plt.show()
