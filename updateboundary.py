import numpy as np
from inv import *


def updateBoundary(X, Y, DZ, Pres, Gamma, pref, fac1, alpha, beta, eps, kappa):
    X = updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa)
    X = updateBoundary_s(X, Y, fac1)
    Pres, Gamma, pmean, dif = correction_Gamma(X, DZ, Pres, Gamma, pref, alpha, kappa)
    X = correction_M(X, DZ, Pres, alpha, kappa)
    X = updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa)
    return X, Pres, Gamma, pmean


def updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa):
    M_int = np.arange(1, X.shape[0], 2)
    # top
    X[M_int, -1] = X[M_int, -3] + 2 * DZ * (1 - alpha) ** kappa / (kappa * alpha)
    # bottom
    t = (Y[::2] ** 2 - X[::2, 1] ** 2) / (1 - X[::2, 1] ** 2)
    t[0] = 0
    t[-1] = 0
    t3 = eps * DZ * (t[:-1] + t[1:]) / (8 * beta)
    X[M_int, 0] = X[M_int, 2] - 2 * DZ / beta * X[M_int, 1] + t3
    return X


def updateBoundary_s(X, Y, fac1):
    X = updateBoundary_s_Iter(X, Y, fac1, 0)
    X = updateBoundary_s_Iter(X, Y, fac1, X.shape[1] - 1)
    return X


def updateBoundary_s_Iter(X, Y, fac1, k):
    Sl_int = np.arange(2, X.shape[1] - 1, 2)
    if k == 0:
        Xold = X[Sl_int, k + 1]
    else:
        Xold = X[Sl_int, k - 1]

    for it in range(4):
        Gkp = (
            fac1[Sl_int] * ((Xold**2 - Y[Sl_int] ** 2) / (1 - Xold**2))
            + X[Sl_int + 1, k]
            - X[Sl_int - 1, k]
        )
        Gkpp = fac1[Sl_int] * 2 * Xold * (1 - Y[Sl_int] ** 2) / (1 - Xold**2) ** 2
        X[Sl_int, k] = Xold - Gkp / Gkpp
        Xold = X[Sl_int, k]
    return X


def correction_Gamma(X, DZ, Pres, Gamma, pref, alpha, kappa):
    M_int = np.arange(1, X.shape[0] - 1, 2)

    Pres[M_int, 2:-2] = (
        alpha * kappa * (X[M_int, 3:-1] - X[M_int, 1:-3]) / (2 * DZ)
    ) ** (1 / kappa)
    pmean = np.sum(Pres[M_int, :] * (X[M_int + 1, :] - X[M_int - 1, :]), axis=0) / 2
    dif = pmean - pref
    Pres[M_int, 1:-1] = Pres[M_int, 1:-1] - dif[:, 1:-1]
    Gamma[M_int, 1:-1] = Pres[M_int, 1:-1] ** (kappa - 1)
    pmean = np.sum(Pres[M_int, :] * (X[M_int + 1, :] - X[M_int - 1, :]), axis=0) / 2
    dif = pmean - pref
    return Pres, Gamma, pmean, dif


def correction_M(X, DZ, Pres, alpha, kappa):
    # Jt = X.shape[0]
    for j in range(1, X.shape[0], 2):
        Mp = Pres[j, 1:-1] ** kappa / (kappa * alpha)
        X[j, 1:-1] = intdde(Mp, X[j, 1], DZ)
    return X
