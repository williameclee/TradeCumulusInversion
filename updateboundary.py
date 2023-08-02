import numpy as np
import inv

def m_boundary(x, scs, dz, tdz, alpha, beta, eps, kappa):
    # top
    mtop = x[::2, -1] = x[::2, -3] + tdz * (1 - alpha) ** kappa / (kappa * alpha)
    # bottom
    t1 = np.zeros_like(x[::2, 0])
    t2 = np.zeros_like(x[::2, 0])

    xsp = x[1, 1] ** 2
    t1[0] = 0
    t2[0] = (xsp - scs[1]) ** 2 / (1 - xsp)
    xsm = x[-2, 1] ** 2
    t1[-1] = (xsm - scs[-2]) ** 2 / (1 - xsm)
    t2[-1] = 0
    xsp = x[3:-1:2, 1] ** 2
    xsm = x[1:-3:2, 1] ** 2
    t1[1:-1] = (xsm - scs[1:-4:2]) ** 2 / (1 - xsm)
    t2[1:-1] = (xsp - scs[4:-2:2]) ** 2 / (1 - xsp)
    term3 = eps * dz * (t1 + t2) / (8 * beta)
    mbottom = x[::2, 2] - tdz / beta * x[::2, 1] + term3

    return mtop, mbottom

def s_boundary(x, fac1, scs):
    for iter in range(4):
        # xold = x[1::2, -1]
        # xs = xold**2
        # gkp = fac1[2:-1:2] * ((xs - scs[2:-1:2]) / (1 - xs)) + x[2::2, -1] - x[:-1:2, -1]
        # gkpp = fac1[2:-1:2] * 2 * xold * (1 - scs[2:-1:2]) / (1 - xs) ** 2
        # stop = xold - gkp / gkpp
        stop = s_hor(x[1::2, -2], x, -1, fac1, scs)
        sbottom = s_hor(x[1::2, 1], x, 0, fac1, scs)
    return stop, sbottom

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
    x[::2, -1], x[::2, 0] = m_boundary(x, scs, dz, tdz, alpha, beta, eps, kappa)
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
    return pres[1:-1:2,:], gamma[1:-1:2,:], avep, dif