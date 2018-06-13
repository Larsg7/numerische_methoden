import numpy as np
import sys
import matplotlib.pyplot as plt


class Integration:
    def __init__(self, N, start, end, f, dim):
        self.f = f
        self.N = N
        self.start = start
        self.end = end
        self.dim = dim
        self.volume = np.abs(self.start - self.end) ** self.dim

    def _gen_random(self) -> np.array:
        return np.random.uniform(self.start, self.end, (self.N, self.dim))

    def calc(self) -> float:
        result = sum((self.f(x)
                      for x in self._gen_random()), 0.0) / self.N

        f_2 = sum((self.f(x) ** 2 for x in self._gen_random()), 0.0) / self.N
        error = self.volume * np.sqrt(np.abs(f_2 - result ** 2) / (self.N - 1))
        return self.volume * result, error


def f(x):
    distance = np.linalg.norm(x)
    return 1 if distance <= 1 else 0


ds = [2, 3, 10]
Ns = [10000, 1000000]
analytical = {
    2: np.pi,
    3: 4/3 * np.pi,
    10: np.pi ** 5 / 120
}

for d in ds:
    for N in Ns:
        i = Integration(N, -1, 1, f, d)
        result, error = i.calc()
        print(f'N = {N:>9}, d = {d:>2}: {result:f} +- {error:1.8f}')

    print(f'Analytical result:     {analytical[d]:f}')
    print()

'''
N =     10000, d =  2: 3.132400 +- 0.01657331
N =   1000000, d =  2: 3.142360 +- 0.00164063
Analytical result:     3.141593

N =     10000, d =  3: 4.215200 +- 0.04006407
N =   1000000, d =  3: 4.185744 +- 0.00400482
Analytical result:     4.188790

N =     10000, d = 10: 2.355200 +- 0.46868777
N =   1000000, d = 10: 2.467840 +- 0.05060461
Analytical result:     2.550164
'''
