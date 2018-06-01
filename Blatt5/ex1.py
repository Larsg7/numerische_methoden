import numpy as np
import sys
import matplotlib.pyplot as plt

n = int(sys.argv[1])
R = 1
a = 2*R/(2*n+1)
matrix_size = (2*n + 1, 2*n + 1)
A = np.zeros(matrix_size)
b = np.zeros(2*n + 1)
b[int(2 * n/2)] = 1/a

def bootstrap_matrix(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            if i == j:
                A[i][j] = -2
            elif (j - 1 == i) or (j + 1 == i):
                A[i][j] = 1
    return 1/a**2 * A


A = bootstrap_matrix(A)
A[0][1] = 0
A[matrix_size[0] - 1][matrix_size[1] - 2] = 0

def lu_decomposition(A: np.array):
    size = A.shape
    L = np.identity(size[0])
    U = np.zeros(size)

    for k in range(size[1]):
        for j in range(k + 1):
            sum = 0
            for l in range(j):
                sum += L[j][l] * U[l][k]
            U[j][k] = A[j][k] - sum
        
        for j in range(k + 1, size[0]):
            sum = 0
            for l in range(k):
                sum += L[j][l] * U[l][k]

            L[j][k] = 1/(U[k][k]) * (A[j][k] - sum)

    return L, U

def calc_L_y(L: np.array, b: np.array):
    size = L.shape[0]
    y = np.zeros(size)
    for j in range(size):
        y[j] = b[j] - sum(L[j][k]*y[k] for k in range(0, j))

    return y

def calc_U_x(U: np.array, y: np.array):
    size = U.shape[0]
    x = np.zeros(size)
    for j in range(size - 1, -1, -1):
        s = sum(U[j][k] * x[k] for k in range(j+1, size))
        x[j] = (y[j] - s) / U[j][j]

    return x


def calc(A, b):
    L, U = lu_decomposition(A)
    y = calc_L_y(L, b)
    return calc_U_x(U, y)

phi = calc(A,b)

x = np.linspace(-R, R, 2*n + 1)

# print(A)
print(phi)


plt.plot(x, phi)
plt.show()