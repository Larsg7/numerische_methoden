from abc import ABC, abstractmethod
import numpy as np


class RungeKutta_4:
    def __init__(self, f, h: float):
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

    def advance_y(self, y: np.array, t: float):
        return y + 1 / 6 * (self.__k1(y, t) + 2 * self.__k2(y, t) + 2 * self.__k3(y, t) + self.__k4(y, t))


class Simulation (ABC):
    @staticmethod
    @abstractmethod
    def force(y, t) -> np.array:
        pass

    def __init__(self, dt):
        self.dt = dt
        self._init_rungekutta()
        self.dts = []
        self.ts = []

    def _init_rungekutta(self):
        self.rungekutta = RungeKutta_4(self._f, self.dt)

    def _f(self, y, t):
        return np.array([y[1], self.force(y[0], t)])

    def _do_step(self, y: np.array, t: float) -> np.array:
        return self.rungekutta.advance_y(y, t)

    def run(self, t_max, y0):
        t = 0
        next_result = self._do_step(y0, t)
        result = []
        while t <= t_max:
            self.dts.append(self.dt)
            result.append(next_result[0])
            next_result = self._do_step(next_result, t)
            t += self.dt
            self.ts.append(t)
        return result


class SimulationAdaptive (Simulation):
    def __init__(self, dt, max_error):
        self.max_error = max_error
        super().__init__(dt)

    def _init_rungekutta(self):
        self.rungekutta = RungeKutta_4(self._f, self.dt)
        self.rungekutta_half = RungeKutta_4(self._f, self.dt / 2)

    def _do_step(self, y_in, t):
        while True:
            self._init_rungekutta()
            y = self.rungekutta.advance_y(y_in, t)
            y_half = self.rungekutta_half.advance_y(y_in, t)
            error = self.est_error(y, y_half)
            self.dt = np.clip(self.dt * 0.9 * (self.max_error / error)
                              ** (1/(4+1)), 0.2 * self.dt, 5 * self.dt) if error != 0 else 1
            if error < self.max_error:
                return y_half

    @staticmethod
    @abstractmethod
    def est_error(y, y_half):
        pass
