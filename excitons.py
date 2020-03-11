print("Calculate Hamiltinian")
h_time = time.time()
mu = np.zeros((1, 3))
mu[0, 0] = 5.09777
mu[0, 1] = 0.07312
mu[0, 2] = 0.04931
A = 219500  # Hartree to cm-1
bohr3 = 0.529177210  # bohr in Å
H = np.zeros((len(cl.molecules), len(cl.molecules)))
for n1 in range(len(cl.molecules)):
    for m1 in range(n1 + 1, len(cl.molecules)):
        v1_o = np.zeros((1, 3))
        v2_o = np.zeros((1, 3))
        v1_o[0, 0] = cl.mass_centers[n1, 0]
        v1_o[0, 1] = cl.mass_centers[n1, 1]
        v1_o[0, 2] = cl.mass_centers[n1, 2]
        v2_o[0, 0] = cl.mass_centers[m1, 0]
        v2_o[0, 1] = cl.mass_centers[m1, 1]
        v2_o[0, 2] = cl.mass_centers[m1, 2]
        dv1 = v1_o - v2_o
        r1 = np.linalg.norm(dv1)
        r1 = r1 / bohr3
        r3 = r1 * r1 * r1
        r5 = r3 * r1 * r1
        J = np.inner(mu, mu) / r3 + (np.inner(mu, dv1) * np.inner(dv1, mu)) / r5
        H[n1, m1] = J * A
        H[m1, n1] = H[n1, m1]
print("   Done: %s" % (time.time() - h_time))
mu_a = np.zeros((len(cl.molecules), 3))
for x in range(len(cl.molecules)):
    mu_a[x, 0] = mu[0, 0]
    mu_a[x, 1] = mu[0, 1]
    mu_a[x, 2] = mu[0, 2]
E, c = np.linalg.eig(H)
bins = 1000
sigma = 100000
EminA = np.min(E)
EmaxA = np.max(E)
print(EmaxA)
Emin = EminA - 0.1 * (EmaxA - EminA) - 3 * sigma
Emax = EmaxA + 0.1 * (EmaxA - EminA) + 3 * sigma
dE = (Emax - Emin) / bins
Ex = np.linspace(Emin, Emax, bins)
Ey = np.zeros(bins)
for n1 in range(len(cl.molecules)):
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
plt.plot(Ex, Ey / np.max(Ey))
# plt.plot(Ex,Cy/np.max(Cy))
plt.plot(Ex, Ny / np.max(Ny))
plt.show()