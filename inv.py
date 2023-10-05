import numpy as np
from scipy.interpolate import interp1d


# def intdde(n, h, u0, f):
#     u = np.zeros(n + 1)
#     # u = np.zeroslike(f)

#     # if method == 1:  # Cubic Interpolation
#     hby24 = h / 24
#     u[0] = u0
#     u[1] = u0 + hby24 * (9.0 * f[0] + 19 * f[1] - 5 * f[2] + f[3])
#     for i in range(2, n):
#         u[i] = u[i - 1] + hby24 * (-f[i - 2] + 13 * (f[i - 1] + f[i]) - f[i + 1])
#     u[n] = u[n - 1] + hby24 * (f[n - 3] - 5 * f[n - 2] + 19 * f[n - 1] + 9 * f[n])
#     # else:  # Trapezoidal Rule
#     #     hby2 = 0.5 * h
#     #     u[0] = u0

#     #     for i in range(1, n + 1):
#     #         u[i] = u[i - 1] + hby2 * (f[i - 1] + f[i])

#     return u

def intdde(Fp, F0,Dx):
    F = np.zeros_like(Fp)
    hby24 = Dx / 24
    F[0] = F0
    F[1] = F0 + hby24 * (9.0 * Fp[0] + 19 * Fp[1] - 5 * Fp[2] + Fp[3])
    for i in range(2, F.shape[0] - 1):
        F[i] = F[i - 1] + hby24 * (-Fp[i - 2] + 13 * (Fp[i - 1] + Fp[i]) - Fp[i + 1])
    F[-1] = F[-2] + hby24 * (Fp[-4] - 5 * Fp[-3] + 19 * Fp[-2] + 9 * Fp[-1])
    return F


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


def pmean(p, x):
    pref = np.sum((p[1:] + p[:-1]) * (x[1:] - x[:-1]) / 2)
    pref = pref / 2
    return pref


# def gtdmsf(B, D, T):
#     n = B.shape[0]
#     for k in range(1, n):
#         B[k,:] = B[k,:] / D[k-1,:]
#         D[k,:] = D[k,:] - B[k,:] * T[k-1,:]
#     return B, D, T


# def gtdmss(B, D, T, X):
#     n = B.shape[0]
#     # Do the forward substitutions
#     for k in range(1, n):
#         X[k,:] = X[k,:] - B[k,:] * X[k-1,:]
#     # Do the backward substitutions
#     X[-1,:] = X[-1,:] / D[-1,:]
#     for k in range(n - 2, -1, -1):
#         X[k,:] = (X[k,:] - T[k,:] * X[k+1,:]) / D[k,:]
#     return X

# def gtm(D,L,U,R,n):
#     A = np.zeros((n, n))
#     d_i, d_j = np.diag_indices_from(A)
#     A[d_i, d_j] = D
#     A[d_i[:-1] + 1, d_j[:-1]] = L[1:]
#     A[d_i[:-1], d_j[:-1] + 1] = U[:-1]

#     X = np.matmul(np.linalg.inv(A),R)
#     return X
