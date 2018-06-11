import numpy as np
import matplotlib.pyplot as plt



def calc_red(A, L, U, j, k, ceil):
    sum = 0
    for l in range(ceil + 1):
        sum += L[j][l] * U[l][k]
    return A[j][k] - sum

def switch(A, k, i_max):
    B = np.copy(A[k])
    C = np.copy(A[i_max])
    A[k] = np.copy(C)
    A[i_max] = np.copy(B)

def pivot(A, L, U, k, b):
    n = len(A)
    l = []
    for i in range(k, n):
        l.append((calc_red(A, L, U, i, k, i-1), i))

    biggest = max(l, key=lambda x: x[0])
    i_max = biggest[1]

    switch(A, k, i_max)
    switch(b, k, i_max)
    switch(L, k, i_max)


def lu_decomposition(A: np.array, b, enable_pivoting = True):
    n =len(A)
    L = np.identity(n)
    U = np.zeros(((n,n)))

    for k in range(n):
        for j in range(k + 1):
            if j == k and enable_pivoting:
                pivot(A, L, U, k, b)
                print('test')
            U[j][k] = calc_red(A, L, U, j, k, j-1)
        
        for j in range(k + 1, n):
            L[j][k] = 1/(U[k][k]) * calc_red(A, L, U, j, k, k-1)

    return L, U

def calc_L_y(L: np.array, b: np.array):
    n = len(L)
    y = np.zeros(n)
    for j in range(n):
        y[j] = b[j] - sum(L[j][k]*y[k] for k in range(0, j))

    return y

def calc_U_x(U: np.array, y: np.array):
    n = len(U)
    x = np.zeros(n)
    for j in range(n - 1, -1, -1):
        s = sum(U[j][k] * x[k] for k in range(j+1, n))
        x[j] = (y[j] - s) / U[j][j]

    return x


def calc(A, b, enable_pivoting = True):
    L, U = lu_decomposition(A, b, enable_pivoting)
    y = calc_L_y(L, b)
    return calc_U_x(U, y)

def ex_i(N):
    A = np.random.uniform(-1, 1, (N, N))
    b = np.random.uniform(-1 ,1, N)
    return np.linalg.norm(A.dot(calc(A, b)) - b), np.linalg.norm(A.dot(calc(A, b, enable_pivoting=False)) - b) 

def ex_ii(n):
    A = np.exp(np.random.uniform(-5, 5, (n, n)))
    b = np.exp(np.random.uniform(-5, 5, n))

    return np.linalg.norm(A.dot(calc(A, b)) - b), np.linalg.norm(A.dot(calc(A, b, enable_pivoting=False)) - b) 

# A = np.array([[1, 1, 1],
#               [7, 3, 0],
#               [2, 4, 3]])

# b = np.array([6, 13, 19])


P = 20
error_i = []
error_i_no_piv = []
error_ii = []
error_ii_no_piv = []
x = []

for i in range(10, P):
    a,b = ex_i(i)
    error_i.append(a)
    error_i_no_piv.append(b)
    a, b = ex_ii(i)
    error_ii.append(a)
    error_ii_no_piv.append(b)
    x.append(i)

plt.plot(x, error_i, label="piv")
plt.show()
print(1)
plt.plot(x, error_i_no_piv, label="no-piv")
plt.show()
print(2)



plt.plot(x, error_ii, label="piv")
plt.show()
print(3)
plt.plot(x, error_ii_no_piv, label="no piv")
plt.show()
print(4)

