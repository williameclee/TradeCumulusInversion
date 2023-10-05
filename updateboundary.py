import numpy as np
import inv


def updateBoundary_M(X, Y, DZ, eps, alpha, beta, kappa):
    M_int = np.arange(1, X.shape[0], 2)
    # top
    X[M_int, -1] = X[M_int, -3] + 2 * DZ * (1 - alpha) ** kappa / (kappa * alpha)
    # bottom
    t = (Y[::2] ** 2 - X[::2, 1] ** 2) / (1 - X[::2, 1] ** 2)
    t[0] = 0
    t[-1] = 0
    t3 = eps * DZ * (t[0:-1] + t[1:]) / (8 * beta)
    X[M_int, 0] = X[M_int, 2] - 2 * DZ / beta * X[M_int, 1] + t3
    return X


def updateBoundary_s(X, Y, fac1):
    for iter in range(4):
        # xold = x[1::2, -1]
        # xs = xold**2
        # gkp = fac1[2:-1:2] * ((xs - scs[2:-1:2]) / (1 - xs)) + x[2::2, -1] - x[:-1:2, -1]
        # gkpp = fac1[2:-1:2] * 2 * xold * (1 - scs[2:-1:2]) / (1 - xs) ** 2
        # stop = xold - gkp / gkpp
        stop = s_hor(X[1::2, -2], X, -1, fac1, Y)
        sbottom = s_hor(X[1::2, 1], X, 0, fac1, Y)
    return stop, sbottom

def updateBoundary_s_Iter(X, Y, fac1, k):
    Sl_int = np.arange(2, X.shape[1]-1, 2)
    if k == 0:
        Xold = X[Sl_int, k+1]
    else:
        Xold = X[Sl_int, k-1]

        for it in range(4):
            Gkp = fac1[Sl_int] * ((Xold ** 2 - Y[Sl_int]**2) / (1 - Xold ** 2)) + X[Sl_int+1, k] - X[Sl_int-1, k]


def s_hor(xold, x, k, fac1, scs):
    for iter in range(4):
        xs = xold**2
        gkp = fac1[2:-1:2] * ((xs - scs[2:-1:2]) / (1 - xs)) + x[2::2, k] - x[:-1:2, k]
        gkpp = fac1[2:-1:2] * 2 * xold * (1 - scs[2:-1:2]) / (1 - xs) ** 2
        xold = xold - gkp / gkpp
    return xold


def m_correction(x, pres, jc, kc, dz, tdz, alpha, beta, eps, kappa, scs):
    # compute a corrected m field
    method = 1
    for j in range(0, 2 * jc - 1, 2):
        ftmp = pres[j, :] ** kappa / (kappa * alpha)
        xtmp = inv.intdde(method, kc, dz, x[j, 0], ftmp)
        x[j, 1:-1] = xtmp
    x[::2, -1], x[::2, 0] = updateBoundary_M(x, scs, dz, tdz, alpha, beta, eps, kappa)
    return x[::2, :]


def gamma_correction(pres, gamma, x, kc, refp, kappa, alpha, tdz):
    # Find gamma for general case
    pres[1:-1:2, :] = (alpha * kappa * (x[::2, 2:] - x[::2, :-2]) / tdz) ** (1 / kappa)
    avep = np.zeros(kc + 1)
    for k in range(kc + 1):
        xp = np.append(x[1:-1:2, k], [1])
        xm = np.append([-1], x[1:-1:2, k])
        avep[k] = np.sum(pres[1:-1:2, k] * (xp - xm)) / 2

    dif = avep - refp
    pres[1:-1:2, :] = pres[1:-1:2, :] - dif[:].T
    for k in range(kc + 1):
        xp = np.append(x[1:-1:2, k], [1])
        xm = np.append([-1], x[1:-1:2, k])
        avep[k] = np.sum(pres[1:-1:2, k] * (xp - xm)) / 2
    gamma[1:-1:2, :] = pres[1:-1:2, :] ** (kappa - 1)
    dif = avep - refp
    return pres[1:-1:2, :], gamma[1:-1:2, :], avep, dif
