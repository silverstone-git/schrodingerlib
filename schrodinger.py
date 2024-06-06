import numpy as np
from sys import argv
import matplotlib.pyplot as plt
#from scipy.constants import hbar

# mass
m = 1
# to ease computation, units are adjusted such that planck's constant = 2*pi
hbar = 1


def scale_kinetic_coeff(el, h):
    #return el*(-1*(hbar**2)/(2*m*h**2))
    return el/h**2


def hamiltonian_gen(n, stencil, h):

    template_matrix = []

    scaled_val = stencil[1]
    for i in range(n):

        cur_row = []
        potential = potential_at(a+(i+1)*h)
        #print("potential being evaluated at: ", a+(i+1)*h)

        # stencil must be refreshed everytime because potential is diff for each row
        #print()
        #print("stencil before potential addition: ", stencil)
        stencil[1] = scaled_val + potential
        #print("stencil after potential addition: ", stencil)
        #print()

        for j in range(i-1):
            cur_row.append(0)

        def step_func(x):
            # centered at 0
            if x < 0:
                return 0
            else:
                return 1

        # putting in the subarray of stencil, ie, 2 elements only in first and last row
        cur_row += stencil[step_func(-i):3-step_func(i-n+1)]

        for j in range(n-i-2):
            cur_row.append(0)

        template_matrix.append(cur_row)

    #for row in template_matrix:
        #print(row)
    res = np.array(template_matrix)
    return res


def potential_at(x):
    return 0




def particle_in_a_box(N, a, b):

    """
    particle in a box using finite element method
    N is the number of wavefunctions in total which are to be considered on the x - axis, from 0th to N-1th
    n is the number of sections in between
    a is x value of first point of the infinite potential well
    b is x value of last point of the infinite potential well
    """

    n = N-1

    h = (b - a)/n

    """
    example -
    if input (N) = 6
    and a = -10,
    b = 10
    points -> 1, 2, 3, 4, 5, 6
    wavefunctions -> psi1, psi2, psi3, psi4, psi5, psi6
    psi1 = psi6 = 0; ignored
    number of sections = 5
    so, step size = x_final - x_init / number of sections
    => h = (10 - (-10)) / 5
         = 20 / 5 = 4
    matrix will only be for psi1, psi2, psi3, psi4, psi5 ie, n = 4
    and thus potentials will be calculated at x = -6, -2, 2, 6

    this example is just for explanation, give N = 10000, a = -50, b = 50 for a good enough approximation
    """

    #print("h is: ", h)

    stencil = [1, -2, 1]
    # scale kinetic operator matrix's coefficient for every element
    stencil = list(map(scale_kinetic_coeff, stencil, [h] * len(stencil)))
    #print("stencil after scaling: ", stencil)

    hamiltonian = hamiltonian_gen(n-1, stencil, h)

    print("hamiltonian is -")
    print(hamiltonian)

    eigenvals, eigenvecs = np.linalg.eig(hamiltonian)

    #for i in eigenvals:
    #    print(i)

    #print("eigenvals are -")
    #print(eigenvals)

    #for i in eigenvecs:
    #    print(i)

    #print(eigenvecs[0])
    xarray = np.array(list([a+h*x for x in range(1, n)]))
    #print(xarray)
    #print(len(xarray))
    #print(len(eigenvecs[0]))
    for eigenvec in eigenvecs:
        yarray = []
        for i in eigenvec:
            yarray.append(i**2)
        # plot of psi squared vs x
        plt.plot(xarray, yarray)
        #plt.plot(xarray, eigenvec)
        plt.show()


    #print("tbd: scaling and adding potential, and then plotting the eigen things")


if __name__ == '__main__':
    a = int(argv[1])
    b = int(argv[2])
    N = int(argv[3])
    particle_in_a_box(N, a, b)
