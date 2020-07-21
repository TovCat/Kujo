import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.special import wofz


class Graph:

    def __init__(self):
        self.data_x = []
        self.data_y = []
        self.data_x_cut = []
        self.data_y_cut = []
        self.peaks = []
        self.peak_threshold = 0.01

    def read(self, path: str):
        contents = []
        try:
            file = open(path, "r")
            contents = file.readlines()
            file.close()
        except OSError:
            print("Could not open the .dat file at: ", path)
            exit(-1)
        for x in contents:
            words = x.split()
            words[0].replace(",", ".")
            words[1].replace(",", ".")
            try:
                self.data_x.append(float(words[0]))
                self.data_y.append(float(words[1]))
            except ValueError:
                print("Erroneous data at: ", path)
                exit(-1)

    def find_peaks(self):
        for i in range(len(self.data_x) - 1):
            dy = self.data_y[i + 1] - self.data_y[i]
            dx = self.data_x[i + 1] - self.data_x[i]
            if abs(dy / dx) <= self.peak_threshold:
                flag = False
                for i1 in range(len(self.peaks)):
                    diff = abs(self.peaks[i1] - i)
                    if diff <= 3:
                        flag = True
                        break
                if not flag:
                    self.peaks.append(i)

    def cut_peak(self, x):
        x_peaks_shifted = []
        for i in range(len(self.peaks)):
            x_peaks_shifted[i] = self.data_x[self.peaks[i]] - x
        peak_of_interest = x_peaks_shifted.index(min(x_peaks_shifted))
        data_x_shifted = []
        for i in range(len(self.data_x)):
            data_x_shifted.append(self.data_x[i] - x)
        min_lim = data_x_shifted[peak_of_interest - 1]
        max_lim = data_x_shifted[peak_of_interest + 1]
        for i in range(len(data_x_shifted)):
            if data_x_shifted[i] > min_lim and data_x_shifted < max_lim:
                self.data_x_cut.append(data_x_shifted[i])
                self.data_y_cut.append(self.data_y[i])

    def fitting(self):

        def voight(x, alpha, gamma):
            sigma = alpha / np.sqrt(2 * np.log(2))
            max_v = np.real(wofz((1j * gamma) / sigma / np.sqrt(2))) / sigma / np.sqrt(2 * np.pi)
            return np.real(wofz((x + 1j * gamma) / sigma / np.sqrt(2))) / sigma / np.sqrt(2 * np.pi) / max_v

        a, b = scipy.optimize.curve_fit(voight, self.data_x_cut, self.data_y_cut)
        return a, b


def generate_plot(E, c, mu, bins, sigma):
    EminA = np.min(E)
    EmaxA = np.max(E)
    if (EmaxA-EminA) / sigma > bins:
        bins = int(np.floor((EmaxA - EminA) / sigma))
    Emin = EminA - 0.1 * (EmaxA - EminA) - 3 * sigma
    Emax = EmaxA + 0.1 * (EmaxA - EminA) + 3 * sigma
    dE = (Emax - Emin) / bins
    Ex = np.linspace(Emin, Emax, bins)
    Ey = np.zeros(bins)
    for n in range(len(E)):
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
    plt.xlabel('eV')
    plt.ylabel('Intensity')
    plt.savefig('plot.png')
