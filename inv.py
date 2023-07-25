import numpy as np
from scipy.interpolate import interp1d


def intdde(method, n, h, u0, f):
    u = np.zeros(n + 1)
    # u = np.zeroslike(f)

    if method == 1:  # Cubic Interpolation
        hby24 = h / 24
        u[0] = u0
        u[1] = u0 + hby24 * (9.0 * f[0] + 19 * f[1] - 5 * f[2] + f[3])
        for i in range(2, n):
            u[i] = u[i - 1] + hby24 * (-f[i - 2] + 13 * (f[i - 1] + f[i]) - f[i + 1])
        u[n] = u[n - 1] + hby24 * (f[n - 3] - 5 * f[n - 2] + 19 * f[n - 1] + 9 * f[n])
    else:  # Trapezoidal Rule
        hby2 = 0.5 * h
        u[0] = u0

        for i in range(1, n + 1):
            u[i] = u[i - 1] + hby2 * (f[i - 1] + f[i])

    return u


def intp(x, sc, ntim):
    # interpolate values to latitude space
    jc2 = jc // 2

    ppos = np.zeros((2 * jl + 1, kc + 1, 21))  # (-jl:jl,0:kc,0:1)
    tpos = np.zeros((2 * jl + 1, kc + 1, 21))  # (-jl:jl,0:kc,0:1)

    plat = np.zeros((2 * jc2 + 1, kc + 1))  # (-jc2:jc2,0:kc)
    pllat = np.rad2deg(np.arcsin(sc[1:-1:2]))
    pllat[0] = -90.0
    pllat[-1] = 90.0

    if ntim == 0:
        plat = np.tile(pllat, (kc + 1, 1)).T
        ppos[:, :, 0] = spx(plat, pllat, lat)
    else:
        plat[1:-1, :] = np.rad2deg(x[2:-1:2, 1:-1])
        plat[0, :] = -90.0
        plat[-1, :] = 90.0
        ppos[:, :, ntim] = spx(plat, pllat, lat)

    for k in range(kc + 1):
        for j in range(-jl, jl + 1):
            tpos[j, k, ntim] = theta[k]

        # Compute theta position
        for k in range(kc + 1):
            for j in range(-jl, jl + 1):
                tpos[j, k, ntim] = theta[k]
    return ppos, tpos


def spx(f, x, xe):
    n, m = f.shape
    ne = xe.shape[0]
    fe = np.zeros((ne, m))

    # Loop through each field column and perform spline interpolation
    for k in range(m):
        spline_fit = interp1d(x, f[:, k], kind="cubic")
        fe[:, k] = spline_fit(xe)

    return fe


def fexp(x):
    k = 1 / 2 * np.exp(2) * np.log(2)
    fac = np.exp(1 / (x - 1))
    fexp_val = np.exp(-k * fac / x)
    return fexp_val


def gtdmsf(m, n, B, D, T, inc, jump):
    inc = int(inc)
    # Compute the factorisations
    l1 = int(0)
    l2 = int(l1 + (m - 1) * jump)
    for k in range(1, n):
        l1 += inc
        l2 += inc
        for l in range(l1, l2 + 1, jump):
            B[l] = B[l] / D[l - inc]
            D[l] = D[l] - B[l] * T[l - inc]
    return B, D, T


def gtdmss(m, n, B, D, T, X, inc, jump):
    # Do the forward substitutions
    l1 = 0
    l2 = l1 + (m - 1) * jump
    for k in range(1, n):
        l1 += inc
        l2 += inc
        for l in range(l1, l2 + 1, jump):
            X[l] = X[l] - B[l] * X[l - inc]

    # Do the backward substitutions
    for l in range(l1, l2 + 1, jump):
        X[l] = X[l] / D[l]
    for k in range(n - 1, 0, -1):
        l1 -= inc
        l2 -= inc
        for l in range(l1, l2 + 1, jump):
            X[l] = (X[l] - T[l] * X[l + inc]) / D[l]
    return X
