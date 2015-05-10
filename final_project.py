from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

def setup_matrix(N, epsilon, delta_t):
    m = 1
    h = 1
    upper = np.eye(N, k=-1)
    lower = np.eye(N, k=1)
    mid = (-2 + (4j*m*epsilon**2)/(h*delta_t))
    return upper+lower+mid


def find_x(phi, delta_t):
    m = 1
    h = 1
    N = len(phi)
    epsilon = 2/(N - 1)
    A = setup_matrix(N, epsilon, delta_t)
    rhs = ((8j*m*epsilon**2)/(h*delta_t))*phi
    return np.linalg.solve(A, rhs)

def create_gaussian_phi(epsilon, k, x0, sigma):
    num = (2*epsilon**-1) + 1
    xx = np.linspace(1, -1, num)
    phi = np.exp(1j*k*xx)*np.exp(-(xx-x0)**2 / (2*sigma**2))
    return phi

def plot_func(func):
    magnitude = np.multiply(np.conj(func), func)
    plt.plot(magnitude)
    plt.show()


def find_nth_state(phi, delta_t, n):
    for i in range(n):
        x = find_x(phi, delta_t)
        phi = x - phi
    return phi
