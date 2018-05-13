import numpy as np
import matplotlib.pyplot as plt
import functools
from rungekutta import Simulation, SimulationAdaptive

h = 0.001
t_max = 20


class KepplerProblem (Simulation):
    @staticmethod
    def force(y, t):
        return - y / (np.linalg.norm(y) ** 3)


class KepplerProblemAdaptive (SimulationAdaptive):
    @staticmethod
    def force(y, t):
        return - y / (np.linalg.norm(y) ** 3)

    @staticmethod
    def est_error(y, y_half) -> float:
        return (np.linalg.norm(y_half[0] - y[0])) / (2 ** 4 - 1)


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


def plot(r0, v0, color='', adaptive=True, max_error=0.0001):
    y = np.array([r0, v0])

    b = get_b(np.linalg.norm(r0))
    c = 1 - np.linalg.norm(r0)

    aphel = get_aphel(r0)
    min_distance_to_aphel = 100

    sim = KepplerProblemAdaptive(
        h, max_error) if adaptive else KepplerProblem(h)

    result = sim.run(t_max, y)
    dts = sim.dts
    plt.subplot(1, 2, 1)
    plt.plot(sim.ts, dts, label='Time step r0={}'.format(r0))
    plt.legend(loc='upper right')

    for r in result:
        dist_to_aphel = np.linalg.norm(aphel - r[0])
        min_distance_to_aphel = min(dist_to_aphel, min_distance_to_aphel)

    print('Min distance to aphel for r0={}: {}'.format(r0, min_distance_to_aphel))

    plt.subplot(1, 2, 2)
    plt.plot([x[0] for x in result], [x[1] for x in result])
    plt.title('r0 = {}, v0 = {}'.format(r0, v0))
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
    time_starts = [3, 6]
    for t in time_starts:
        index = int(t/h)
        length = int(time / h)
        data = result[index:index + length]
        area = functools.reduce(lambda acc, p: (
            acc[0] + area_centric_triangle(acc[1], p), p), data, (0, result[index]))[0]
        print('The area for time={} and starting time {} is: {}'.format(time, t, area))


def main():
    init_cond = [
        0.9,
        0.1
    ]
    colors = [
        'g',
        'b',
        'y',
        'm'
    ]
    for i, r in enumerate(init_cond):
        v = np.array([0, v0(r)])
        r = np.array([r, 0])

        result = plot(r, v, color=colors[i])

        # verify_second_law(result)
        print()
    plt.title('dt={}, t_max={}'.format(h, t_max))
    plt.show()


main()
