import numpy as np
import matplotlib.pyplot as plt
from eigenvalue import jacobi


def bootstrap_K(N):
    matrix_size = (N ** 2, N ** 2)
    A = np.zeros(matrix_size)

    for i in range(len(A)):
        for j in range(len(A[i])):
            if i == j:
                A[i][j] = 4.
            elif ((j - 1) % (N ** 2) == i) or ((j + 1) % (N ** 2) == i):
                A[i][j] = -1.
            elif ((j - (N)) % (N ** 2) == i) or ((j + (N)) % (N ** 2) == i):
                A[i][j] = -1.

    return A


def bootstrap_K_free(N):
    matrix_size = (N ** 2, N ** 2)
    A = np.zeros(matrix_size)

    for i in range(len(A)):
        for j in range(len(A[i])):
            if i == j:
                a = 4.
                if np.abs(i % N - (N - 1)) % (N - 1) == 0:
                    a -= 1
                if np.abs(int(i / N) - (N - 1)) % (N - 1) == 0:
                    a -= 1
                A[i][j] = a
            elif ((j - 1) == i) or ((j + 1) == i):
                A[i][j] = -1.
            elif ((j - (N)) == i) or ((j + (N)) == i):
                A[i][j] = -1.

    return A


# print([[4 - 2 * np.cos(2 * np.pi * i / 4) - 2 * np.cos(2 * np.pi * j / 4)
#         for i in range(4)] for j in range(4)])

def initial_condition(vs, N):
    V_ = np.linalg.inv(vs)
    x0 = np.zeros(N * N)
    x0[0] = 1
    a = V_.dot(x0)
    return a


# def x(ws, vs, t, j, a):
#     cycl = np.array([])
#     for i in range(len(vs)):
#         cycl = np.append(cycl, (vs[i][j]*np.cos(ws[i] * t) * a[i]), axis=0)

#     return cycl


def x(ws, vs, time, j, a):
    def xt(t):
        return np.array([vs[j][mode] * np.cos(ws[mode] * t) for mode in range(len(ws))]).dot(a)
    return [xt(t) for t in time]


def iv():
    N = 4
    K = bootstrap_K(N)
    print(K)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]
    print('Non-free')
    print(ws)
    print()
    print(vs)
    a = initial_condition(vs, N)
    b = np.array([vs[i][1] for i in range(len(vs))])
    print(b.dot(a))
    return ws, vs


def v():
    print('Free')
    N = 4
    K = bootstrap_K_free(N)
    print(K)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]
    print(ws)
    print()
    print(vs)

    resultBounded = iv()
    ws_ = resultBounded[0]
    vs_ = resultBounded[1]
    a = initial_condition(vs, N)
    a_ = initial_condition(vs_, N)

    j = 15
    time = np.arange(0, 30, 0.1)
    plt.plot(time, x(ws, vs, time, j, a), label="Periodic j = 15")
    plt.plot(time, x(ws_, vs_, time, j, a_), label="Free j = 15")

    plt.legend()
    plt.show()


def vi():
    time = np.arange(0, 30, 0.1)
    J = [0, 48]
    N = 7

    K = bootstrap_K_free(N)
    print(K)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]

    a = initial_condition(vs, N)

    for j in J:
        y = x(ws, vs, time, j, a)
        plt.plot(time, y, label=f'j = {j}')
    plt.legend()
    plt.show()


# iv()
# v()
vi()
