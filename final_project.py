from __future__ import division
import numpy as np
from matplotlib import pyplot as plt


def setup_matrix(N, epsilon, delta_t):
    # Sets up matrix for equation 6.
    m = 1
    h = 1
    upper = np.eye(N, k=-1)
    lower = np.eye(N, k=1)
    mid = (-2 + (4j*m*epsilon**2)/(h*delta_t))*np.eye(N)
    matrix = upper+lower+mid
    return matrix


def find_x(phi, delta_t):
    # Solves equation 6 for X, given phi and delta_t.
    m = 1
    h = 1
    N = len(phi)
    epsilon = 2/(N - 1)
    A = setup_matrix(N, epsilon, delta_t)
    rhs = ((8j*m*epsilon**2)/(h*delta_t))*phi
    return np.linalg.solve(A, rhs)


def create_gaussian_phi(epsilon, k, x0, sigma):
    # Sets up a gaussian using linspace.
    num = (2*epsilon**-1) + 1
    xx = np.linspace(-1, 1, num)
    #phi = np.exp(1j*k*xx)*np.exp(-(xx-x0)**2 / (2*sigma**2))
    #phi = np.exp(-(xx-x0)**2 / (2*sigma**2))
    phi = (np.cos(k*xx) + 1j*np.sin(k*xx))*(np.exp(-(xx-x0)**2 / (2*sigma**2)))
    return phi


def plot_magnitude(func):
    # Plots magnitude of function.
    magnitude = np.multiply(np.conj(func), func)
    plt.plot(magnitude)
    plt.show()


def plot_real_and_imag(func):
    # Plots real and imaginary aspects of function.
    f, subplots = plt.subplots(2, sharex=True)
    subplots[0].plot(np.real(func))
    subplots[0].set_title("Real")
    subplots[1].plot(np.imag(func))
    subplots[1].set_title("Imaginary")
    f.show()


def find_nth_state(phi, delta_t, n):
    # Iterates equation 6 n times.
    for i in range(n):
        x = find_x(phi, delta_t)
        phi = x - phi
    return phi
