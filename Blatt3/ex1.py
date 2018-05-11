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


def v0(r0):
    if r0 == 1:
        return 1
    return np.sqrt((-3*r0+r0**2+2)/(r0 - r0**2))


def get_aphel(r0):
    return np.array([r0[0]-2, r0[1]])


def get_b(r0):
    a = 1
    c = a - r0
    return np.sqrt(a**2 - c**2)


def force(r):
    return - r / (np.linalg.norm(r) ** 3)


def f(y, t):
    return np.array([y[1], force(y[0])])


def plot(r0, v0, color=''):
    y = np.array([r0, v0])

    R = RungeKutta_4(f, h)

    b = get_b(np.linalg.norm(r0))
    c = 1 - np.linalg.norm(r0)

    result = []
    aphel = get_aphel(r0)
    min_distance_to_aphel = 100

    for i in range(int(t_max / h)):
        t = i * h
        result.append(y[0])
        y = R.advance_y(y, t)
        dist_to_aphel = np.linalg.norm(aphel - y[0])
        min_distance_to_aphel = min(dist_to_aphel, min_distance_to_aphel)

    print('Min distance to aphel for r0={}: {}'.format(r0, min_distance_to_aphel))

    plt.plot([x[0] for x in result], [x[1] for x in result])
    # plt.title('r0 = {}, v0 = {}'.format(r0, v0))
    plt.plot([0], [0], 'ro')
    plt.plot(r0[0], r0[1], color+'o', label='Start r0={}'.format(r0))
    plt.plot(r0[0]-2, r0[1], color+'o', label='Aphel r0={}'.format(r0))
    plt.plot(-c, b, color+'o', label='Maxima for minor axis')
    plt.plot(-c, -b, color+'o')
    plt.xlabel('x/a')
    plt.ylabel('y/a')
    plt.legend(loc='upper right')
    return result


def verify_second_law(result):
    def area_centric_triangle(a, b):
        c = np.array([0, 0])
        return 0.5 * np.linalg.norm(np.cross(b-a, c-a))

    time = 1
    time_starts = [6, 11, 18]
    for t in time_starts:
        index = int(t/h)
        length = int(time / h)
        data = result[index:index + length]
        area = functools.reduce(lambda acc, p: (
            acc[0] + area_centric_triangle(acc[1], p), p), data, (0, result[index]))[0]
        print('The area for time={} and starting time {} is: {}'.format(time, t, area))


def main():
    init_cond = [
        1,
        0.5,
        0.3,
        0.1
    ]
    colors = [
        'g',
        'b',
        'y',
        'b'
    ]
    for i, r in enumerate(init_cond):
        v = np.array([0, v0(r)])
        r = np.array([r, 0])

        result = plot(r, v, color=colors[i])

        verify_second_law(result)
        print()
    plt.title('dt={}, t_max={}'.format(h, t_max))
    plt.show()


main()
