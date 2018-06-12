import numpy as np
import sys
import matplotlib.pyplot as plt
from time import sleep

def one_d():
    n = 100
    R = 1
    a = 2*R/(2*n+1)
    matrix_size = (2*n + 1, 2*n + 1)
    A = np.zeros(matrix_size)
    b = np.zeros(2*n + 1)
    b[int(2 * n/2)] = 1/a
    for i in range(len(A)):
        for j in range(len(A[i])):
            if i == j:
                A[i][j] = -2
            elif (j - 1 == i) or (j + 1 == i):
                A[i][j] = 1
    A = 1/a**2 * A

    return CG(A, b, 0.01).run()

class CG:
  def __init__(self, A: np.array, b: np.array, max_error):
    self.A = A
    self.b = b
    self.max_error = max_error
    self.x = np.zeros(b.shape)
    self.r = b - self.A.dot(self.x)
    self.p = self.r

  def _error(self) -> float:
    return np.linalg.norm(self.b - self.A.dot(self.x))

  def _alpha(self) -> float:
    return (self.r.dot(self.r)) / (self.p.dot(self.A).dot(self.p))

  def _min_alpha(self) -> float:
    return (self.b.dot(self.p) - 0.5 * (self.x.dot(self.A).dot(self.p) + self.p.dot(self.A).dot(self.x))) / (self.p.dot(self.A).dot(self.p))

  def _beta(self, r_n: np.array) -> float:
    return (self.r.dot(self.r)) / (r_n.dot(r_n))

  def run(self):
    while self._error() > self.max_error:
      alpha = self._alpha()
      self.x = self.x + alpha * self.p
      r_n = np.copy(self.r)
      self.r = self.r - alpha * self.A.dot(self.p)
      beta = self._beta(r_n)
      self.p = self.r + beta * self.p
      print(self.x, self._error())
    
    return self.x

A = np.array([[1,2], [2,5]])

b = np.array([5,12])

cg = CG(A, b, 0.01)
print(cg.run() == np.array([1,2]))
print()

one_d()

