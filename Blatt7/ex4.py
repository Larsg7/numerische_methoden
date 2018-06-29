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
            elif ((j - (N + 1)) % (N ** 2) == i) or ((j + (N + 1)) % (N ** 2) == i):
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
            elif ((j - (N + 1)) == i) or ((j + (N + 1)) == i):
                A[i][j] = -1.

    return A


# print([[4 - 2 * np.cos(2 * np.pi * i / 4) - 2 * np.cos(2 * np.pi * j / 4)
#         for i in range(4)] for j in range(4)])

def x(ws, vs, time, j):
    return [sum(vs[mode][j] * np.cos(ws[mode] * t) for mode in range(len(ws))) for t in time]


def iv():
    K = bootstrap_K(4)
    print(K)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]
    print('Non-free')
    print(ws)
    print()
    print(vs)
    return ws, vs


def v():
    print('Free')
    K = bootstrap_K_free(4)
    # print(K)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]
    print(ws)
    print()
    print(vs)

    resultBounded = iv()
    ws_ = resultBounded[0]
    vs_ = resultBounded[1]

    time = np.arange(0, 30, 0.1)
    plt.plot(time, x(ws, vs, time, 15), label="Periodic j = 15")
    plt.plot(time, x(ws_, vs_, time, 15), label="Free j = 15")

    plt.legend()
    plt.show()


def vi():
    time = np.arange(0, 30, 0.1)
    J = [0, 6, 24, 48]
    N = 7

    K = bootstrap_K_free(N)
    result = jacobi(K)
    ws = result[0]
    vs = result[1]

    for j in J:
        y = x(ws, vs, time, j)
        plt.plot(time, y, label=f'j = {j}')
    plt.legend()
    plt.show()


v()
vi()
