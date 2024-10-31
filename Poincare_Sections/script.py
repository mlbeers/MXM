# Importing libraries for data processing, mathematical calculations, and surface dynamics
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import os
import pwlf
from surface_dynamics.all import *
from Library import *
from Library import Section
import math
from time import time
import copy
from scipy import integrate
import sympy as sym
from sympy import Symbol, solve, lambdify
import traceback
import dill
import sys
import unittest
from surface_dynamics.all import Origami
from utils import load_arrays_from_file  # testing
from sage.all import matrix  # testing

# Define number of squares and index, input through command line
n_squares = int(sys.argv[1])  # Number of squares for STS
index = int(sys.argv[2])  # Index to start at

# Generate and select a permutation list based on input
permutations = perms_list(n_squares)
perm = permutations[index]

# Load vectors (saddle connections) from precomputed file
vec_file = "vecs" + str(n_squares) + "-" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))
print(f'number of vecs: {len(vecs0)}')

# Generate alpha matrices, transformations, and eigenvectors for analysis
a, c, e, g = [], [], [], []
for num in range(50):
    try:
        gs = generators(perm, vecs0)
        alphas, Cs, C_invs, eigs, Ms, gens, eigenvecs = poincare_details(perm, vecs0, gs)
    except:
        continue
    a.append(alphas)
    c.append(Cs)
    e.append(eigenvecs)
    g.append(generators)
print(f'length of alphas: {len(a)}')

# Save generated data to a file for future use
data = [a, c, e, g]
with open(os.path.join("results", f"{n_squares} - {index}", "setup.dill"), 'wb') as f:
    dill.dump(data, f)

# Define step size and initialize list for covolume values
dx = 0.0005
covolume_list = []

# Loop through each cusp to analyze sections and store results
for j in range(len(a[0])):
    improved = True if j == 0 else False  # Use initial improvement for first section only

    for i in range(len(a)):
        # Get section properties for given parameters
        vecs, x_vals, m0, m1, x0, y0, dx_y, z = setup(
            a[i][j], c[i][j], e[i][j], vecs0, dx, improved)
        print("i = " + str(i), "j = " + str(j))

        # Skip if section value too small
        if float(z) <= float(1/50000):
            print("too small")
            continue

        # Create DataFrame of winning vector at points in the section
        df = winners(vecs, x_vals, m0, m1, y0, dx, dx_y)

        # Plot and save Poincare section
        try:
            plot(df, vecs, c[i][j], j, n_squares, index, test=False)
        except Exception as error:
            print(error)
            continue
        df.to_csv(os.path.join("results", f"{n_squares} - {index}", "df - " + str(j)), index=False)

        # Generate section objects and boundaries for subsections
        sec_list = sec_setup(df, dx_y)
        secs = sec_comp(sec_list, dx)
        with open(os.path.join("results", f"{n_squares} - {index}", "secs - " + str(j) + ".dill"), 'wb') as f:
            dill.dump(secs, f)

        # Calculate and plot PDF for each cusp
        times = time_comp(secs) if improved else [1]
        pdf(list(df["time"]), times, dx * 2, n_squares, index, j)

        # Append covolume calculations to the list
        covolume_list.append(covolume(secs))
        print("section done")
        break

# Write covolume data to a results file
with open(os.path.join("results", f"{n_squares} - {index}", "covolume.txt"), "w") as file:
    total = sum(covolume_list)
    for item in covolume_list:
        file.write(str(item) + "\n")
    file.write("\n" + str(total))