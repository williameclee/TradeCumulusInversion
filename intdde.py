import numpy as np

def intdde(method, n, h, u0, f):
    u = np.zeros(n + 1)
    print(n)

    if method == 1:  # Cubic Interpolation
        hby24 = h / 24.0
        u[0] = u0
        u[1] = u0 + hby24 * (9.0 * f[0] + 19.0 * f[1] - 5.0 * f[2] + f[3])

        for i in range(1, n+1):
            print(i)
            u[i] = u[i - 1] + hby24 * (-f[i - 2] + 13.0 * (f[i - 1] + f[i]) - f[i + 1])
        print(n)
        u[n] = u[n - 1] + hby24 * (f[n - 3] - 5.0 * f[n - 2] + 19.0 * f[n - 1] + 9.0 * f[n])
    else:  # Trapezoidal Rule
        hby2 = 0.5 * h
        u[0] = u0

        for i in range(1, n + 1):
            u[i] = u[i - 1] + hby2 * (f[i - 1] + f[i])

    return u