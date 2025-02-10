import numpy as np

def lennard_jones_potential(r, epsilon, sigma):
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

def lennard_jones_force(r, epsilon, sigma):
    return 24 * epsilon * (2 * (sigma/r)**12 - (sigma/r)**6) / r