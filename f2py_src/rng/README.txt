F2PY Random Number Module


Uses Fortran RANDOM_NUMBER

These may not be the most efficient ways to generate certain distributions, but I had a few hours to write this so...


Example useage in Python

from f_rng import rng

rng.r_float(r)                # single float from 0.0 to 1.0
rng.r_float_1d(r)             # array of floats from 0.0 to 1.0

rng.r_float_ab                # single float from a to b
rng.r_float_ab_1d             # array of floats from a to b


Gaussian:

rng.r_gaussian(x_min, x_max, mu, sigma, r)
rng.r_gaussian_1d(x_min, x_max, mu, sigma, r)


Maxwell-Boltzmann

rng.r_maxwell(x_min, x_max, a, r)
rng.r_maxwell_1d(x_min, x_max, a, r)





