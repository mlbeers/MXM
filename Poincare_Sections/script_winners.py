from multiprocessing import Pool
from sage.all import matrix  # testing
from sage.all import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import os
import pwlf
from surface_dynamics.all import *
import math
from time import time
import copy
from scipy import integrate
import traceback
import dill
import sys
import unittest
from surface_dynamics.all import Origami
from utils import load_arrays_from_file  # testing
from fractions import Fraction as frac
import sympy as sym
from sympy import Symbol
from sympy import solve, lambdify, Eq
from sympy import Rational, sqrt
from Library import *
from Library import Section

t0 = time()

# number of squares for STS
n_squares = int(sys.argv[1])
# index to start at
index = int(sys.argv[2])

dx = float(sys.argv[3])

os.makedirs(os.path.join("results", f"{n_squares}_{index}"), exist_ok=True)  # Ensure directory exists

permutations = perms_list(n_squares)

perm = permutations[index]
plot = perm.plot()
plot.save(os.path.join("results", f"{n_squares}_{index}", "permutation.png"))

# get a list of saddle connections
vec_file = "vecs" + str(n_squares) + "_" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))

print(f'number of vecs: {len(vecs0)}')

# generate a list of alpha, c matrices, and eigenvectors for each cusp of the STS to experiment with to find "nice" sections for our poincare sections
a = []
c = []
e = []
g = []
for num in range(10):
    try:
        gs = generators(perm, vecs0)
        alphas, Cs, Ss, eigs, Ms, gens, eigenvecs = poincare_setup(perm, vecs0, gs)
        print(alphas)
    except Exception as ex:
        print(ex)
        continue
    a.append(alphas)
    c.append(Cs)
    e.append(eigenvecs)
    g.append(gs)
print(f'length of alphas: {len(a)}')

# write these values to a file
data = [a, c, e, g]
with open(os.path.join("results", f"{n_squares}_{index}", "setup.dill"), 'wb') as f:
    dill.dump(data, f)

# get the number of cores to be used in computations and the number of loops needed to complete each cusp
num_pools, num_loops = pool_num(len(a[0]))

args = [(vecs0, a, c, e, dx, dy, dz, j, n_squares, index, estimated) for j in range(len(a[0]))]

# parallelize the run_script function
for i in range(num_loops):
    if i == num_loops - 1:  # Last batch might have fewer tasks
        with Pool(len(args[i * num_pools:])) as p:
            p.starmap(compute_poincare_sections, args[i * num_pools:])
    else:
        with Pool(num_pools) as p:
            p.starmap(compute_poincare_sections, args[i * num_pools : (i + 1) * num_pools])

t1 = time()
print(f"winners done: {(t1-t0)/60**2}\n")
print("-----------------------------------------------------------------------------------")