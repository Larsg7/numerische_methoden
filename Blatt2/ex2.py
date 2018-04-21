import numpy as np
import matplotlib.pyplot as plt

h = 0.001
t_max = 10
phi_start = 1
l = 1
m = 1
g = 9.81

class RungeKutta_4:
    def __init__(self, f, h):
        self.f = f
        self.h = h

    def __k1(self, y, t):
        return self.f(y, t) * self.h

    def __k2(self, y, t):
        return self.f(y + 0.5 * self.__k1(y, t), t + 0.5 * self.h) * self.h

    def __k3(self, y, t):
        return self.f(y + 0.5 * self.__k2(y, t), t + 0.5 * self.h) * self.h

    def __k4(self, y, t):
        return self.f(y + self.__k3(y, t), t + self.h) * self.h

    def advance_y(self, y, t):
        return y + 1 / 6 * (self.__k1(y, t) + 2 * self.__k2(y, t) + 2 * self.__k3(y, t) + self.__k4(y, t))
    

def main():
    def force(phi):
        return - m * g * phi

    def f(y, t):
        return np.array([y[1], force(y[0])])

    def analytical(t):
        return phi_start * np.cos(np.sqrt(g / l) * t)

    y = np.array([phi_start, 0])

    R = RungeKutta_4(f, h)

    result = []
    result_analytical = []
    time = []

    for i in range(int(t_max / h)):
        t = i * h
        time.append(t)
        result.append(y[0])
        result_analytical.append(analytical(t))
        y = R.advance_y(y, t)

    plt.plot(time, result)
    plt.plot(time, result_analytical)
    plt.show()

main()
