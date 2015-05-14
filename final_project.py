from __future__ import division
import numpy as np
from matplotlib import pyplot as plt


def setup_matrix(N, epsilon, delta_t, potential=None):
    # Sets up matrix for equation 6.
    m = 1
    h = 1
    upper = np.eye(N, k=-1)
    lower = np.eye(N, k=1)
    mid = (-2 + (4j*m*epsilon**2)/(h*delta_t))*np.eye(N)
    matrix = upper+lower+mid
    if potential is not None:
        matrix -= ((2*m*epsilon**2)/h**2)*np.diag(potential)
    return matrix


def find_x(phi, delta_t, potential=None):
    # Solves equation 6 for X, given phi and delta_t.
    m = 1
    h = 1
    N = len(phi)
    epsilon = 2/(N - 1)
    A = setup_matrix(N, epsilon, delta_t, potential)
    print A[-1][-1]
    rhs = ((8j*m*epsilon**2)/(h*delta_t))*phi
    return np.linalg.solve(A, rhs)


def find_nth_state(phi, delta_t, n, potential=None):
    # Iterates equation 6 n times.
    for i in range(n):
        x = find_x(phi, delta_t, potential)
        phi = x - phi
    return phi


def create_gaussian_phi(epsilon, k, x0, sigma):
    # Sets up a gaussian using linspace.
    num = (2*epsilon**-1) + 1
    xx = np.linspace(-1, 1, num)
    phi = (np.cos(k*xx) + 1j*np.sin(k*xx))*(np.exp(-(xx-x0)**2 / (2*sigma**2)))
    return phi


def create_step_function(epsilon, h):
    num = (2*epsilon**-1) + 1
    xx = np.linspace(-1, 1, num)
    v =  np.zeros(len(xx))
    v[np.floor(len(xx)/2):] = h
    return v


def plot_phi(func):
    # Plots real and imaginary aspects of function.
    f, subplots = plt.subplots(3, sharex=True)
    xx = np.linspace(-1, 1, len(func))

    magnitude = np.multiply(np.conj(func), func)
    subplots[0].plot(xx, magnitude)
    subplots[0].set_title("Magnitude")
    subplots[0].set_ylabel("Phi*Phi")
    subplots[0].set_xlabel("X")
    subplots[1].plot(xx,np.real(func))
    subplots[1].set_title("Real")
    subplots[1].set_ylabel("Re(Phi)")
    subplots[1].set_xlabel("X")
    subplots[2].plot(xx, np.imag(func))
    subplots[2].set_title("Imaginary")
    subplots[2].set_ylabel("Im(Phi)")
    subplots[2].set_xlabel("X")
    f.show()
