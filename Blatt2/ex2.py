import numpy as np
import matplotlib.pyplot as plt

h = 0.0001
t_max = 5
phi_start = 1
l = 1
m = 1
g = 9.81

phi_starts = [0.01, np.pi / 4, np.pi / 2, 3 * np.pi / 4]

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
    

def get_wave_period(time, result):
    downwards = True
    start_time = 0
    cross_zero = 0
    for index, point in enumerate(result):
        if point < 0 and downwards:
            cross_zero += 1
            downwards = False
        if point > 0 and not downwards:
            downwards = True
        if cross_zero == 1 and start_time == 0:
            start_time = time[index]
        elif cross_zero == 2:
            return time[index] - start_time
        

def force(phi):
        return - m * g * np.sin(phi)

def f(y, t):
    return np.array([y[1], force(y[0])])

def analytical(t, phi):
    return phi * np.cos(np.sqrt(g / l) * t)

def plot(phi):
    y = np.array([phi, 0])

    R = RungeKutta_4(f, h)

    result = []
    result_analytical = []
    time = []

    for i in range(int(t_max / h)):
        t = i * h
        time.append(t)
        result.append(y[0])
        a = analytical(t, phi)
        result_analytical.append(a)
        y = R.advance_y(y, t)

    print('Wave period numeric ({}): {}s'.format(phi, get_wave_period(time, result)))
    print('Wave period analytical ({}): {}s'.format(phi, get_wave_period(time, result_analytical)))

    p1 = plt.plot(time, result)
    p2 = plt.plot(time, result_analytical)
    plt.title('Phi_0 = {}'.format(phi))
    plt.xlabel('time/s')
    plt.ylabel('phi')
    plt.legend((p1[0], p2[0]), ('Numeric result', 'Analytical result'))
    plt.show()

def main():
    for i in phi_starts:
        plot(i)

main()
