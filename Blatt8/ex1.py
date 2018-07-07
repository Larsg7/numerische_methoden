import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

data = """
+1.00000000000e+00 +6.33713107215e-02 3.78369546667e-05
+2.00000000000e+00 +1.44577642493e-01 1.20865422902e-04
+3.00000000000e+00 +2.51824506692e-01 2.56618760202e-04
+4.00000000000e+00 +3.31545011750e-01 4.41171477669e-04
+5.00000000000e+00 +3.98582892282e-01 6.69412637547e-04
+6.00000000000e+00 +4.59482949153e-01 1.05880393429e-03
+7.00000000000e+00 +5.15538066799e-01 1.67403644680e-03
+8.00000000000e+00 +5.68026108073e-01 2.56413283911e-03
+9.00000000000e+00 +6.17446547869e-01 3.82281996304e-03
+1.00000000000e+01 +6.65424320744e-01 5.27638975817e-03
""".strip().split('\n')[2:]

data = [[float(a) for a in d.split(' ')] for d in data]
r = [d[0] for d in data]
V = [d[1] for d in data]


def xi_squared(data, g):
    return sum(((g(r) - x) / o) ** 2 for (r, x, o) in data)


def g(V0, alpha, sigma):
    def f(x):
        return V0 + alpha / x + sigma * x
    return f


def to_minimize(x):
    return xi_squared(data, g(x[0], x[1], x[2]))


def do_minimize():
    res = minimize(to_minimize, np.array([1, 1, 1]))
    print(res.message)
    print(res.x)
    params = res.x
    print(f"Final x^2 / d.o.f: {to_minimize(params) / (len(data) - 4)}")

    return params


def b():
    params = do_minimize()

    plt.plot(r, V, 'o', label="Data points")
    plt.xlabel("r")
    plt.ylabel("V(r)")
    plot_r = np.arange(r[0], 20, 0.01)
    plt.plot(plot_r, [g(params[0], params[1], params[2])(x)
                      for x in plot_r], label="Fit (chi squared)")
    plt.legend()
    plt.show()


def lagrange(data, r):
    N = len(data)

    def l(x, j):
        res = 1
        for k in range(N):
            if k == j:
                continue
            res *= (x - data[k][0]) / (data[j][0] - data[k][0])
        return res

    return sum(data[j][1] * l(r, j) for j in range(N))


def iii():
    params = do_minimize()
    plt.plot(r, V, 'o', label="Data points")
    plt.xlabel("r")
    plt.ylabel("V(r)")
    plot_r = np.arange(2, 12, 0.01)
    plt.plot(plot_r, [lagrange(data, x)
                      for x in plot_r], label="Fit (interpolated)")
    plt.plot(plot_r, [g(params[0], params[1], params[2])(x)
                      for x in plot_r], label="Fit (chi squared)")
    plt.legend()
    plt.show()


# b()
iii()
