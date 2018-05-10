import numpy as np
import matplotlib.pyplot as plt
import functools

h = 0.001
t_max = 20


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


def force(r):
    return - r / (np.linalg.norm(r) ** 3)


def f(y, t):
    return np.array([y[1], force(y[0])])


def plot(r0, v0):
    y = np.array([r0, v0])

    R = RungeKutta_4(f, h)

    result = []

    for i in range(int(t_max / h)):
        t = i * h
        result.append(y[0])
        y = R.advance_y(y, t)

    plt.plot([x[0] for x in result], [x[1] for x in result])
    plt.title('r0 = {}, v0 = {}'.format(r0, v0))
    plt.plot([0], [0], 'ro', label='F')
    plt.plot(r0[0], r0[1], 'go', label='Start')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    return result


def verify_second_law(result):
    def area_centric_triangle(a, b):
        c = np.array([0, 0])
        return 0.5 * np.linalg.norm(np.cross(b-a, c-a))

    time = 1
    time_starts = [6, 11]
    for t in time_starts:
        index = int(t/h)
        length = int(time / h)
        data = result[index:index + length]
        area = functools.reduce(lambda acc, p: (
            acc[0] + area_centric_triangle(acc[1], p), p), data, (0, result[index]))[0]
        print('The area for time={} and starting time {} is: {}'.format(time, t, area))


def main():
    for (r, v) in [(np.array([1, 0]), np.array([0, 1]))]:
        result = plot(r, v)

        verify_second_law(result)
    plt.show()


main()
