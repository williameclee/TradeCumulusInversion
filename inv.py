import numpy as np
from scipy.interpolate import interp1d


def intdde(Fp, F0, Dx):
    F = np.zeros_like(Fp)
    hby24 = Dx / 24
    F[0] = F0
    F[1] = F0 + hby24 * (9 * Fp[0] + 19 * Fp[1] - 5 * Fp[2] + Fp[3])
    for i in range(2, F.shape[0] - 1):
        F[i] = F[i - 1] + hby24 * (-Fp[i - 2] + 13 * (Fp[i - 1] + Fp[i]) - Fp[i + 1])
    F[-1] = F[-2] + hby24 * (Fp[-4] - 5 * Fp[-3] + 19 * Fp[-2] + 9 * Fp[-1])
    return F


def fexp(X):
    gamma = 1 / 2 * np.exp(2) * np.log(2)
    Y = np.exp(-gamma * np.exp(1 / (X - 1)) / X)
    return Y


def pmean(p, x):
    pref = np.sum((p[1:] + p[:-1]) * (x[1:] - x[:-1]) / 2)
    pref = pref / 2
    return pref
