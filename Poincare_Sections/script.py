import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import numpy as np
from matplotlib import pyplot as plt
import os
import pwlf
import os
from surface_dynamics.all import *
from Library import *
from Library import Section
import math
from time import time
import copy
from scipy import integrate
import sympy as sym
from sympy import Symbol
from sympy import solve, lambdify
import traceback
import dill
import sys

# number of squares for STS
n_squares = int(sys.argv[1])
# index to start at
index = int(sys.argv[2])

permutations = perms_list(n_squares)

print("\n\n INDEX - " + str(index) + "\n\n")

perm = permutations[index]

# get a list of saddle connections
vec_file = "vecs" + str(n_squares) + "-" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))

print(len(vecs0))

# generate a list of alpha, c matrices, and eigenvectors for each cusp of the STS to experiment with to find "nice" sections for our poincare sections
a = []
c = []
e = []
g = []
for num in range(200):
    try:
        generators = generators(perm, vecs0)
        alphas, Cs, C_invs, eigs, Ms, generators, eigenvecs = poincare_details(perm, vecs0, generators)
    except:
        continue
    a.append(alphas)
    c.append(Cs)
    e.append(eigenvecs)
    g.append(generators)
print(len(a))

# write these values to a file
data = [a, c, e, g]
with open(os.path.join("results", f"{n_squares} - {index}", "setup.dill"), 'wb') as f:
    dill.dump(data, f)

dx = 0.0005
covolume_list = []

# go thorugh all cusps
for j in range(len(a[0])):

    # change back to False in future
    improved = True
    if j == 0:
        improved = True

    for i in range(len(a)):

        # get dimensions of section
        vecs, x_vals, m0, m1, x0, y0, dx_y, z = setup(
            a[i][j], c[i][j], e[i][j], vecs0, dx, improved)
        print(z)
        print("i = " + str(i), "j = " + str(j))

        if float(z) <= float(1/50000):
            print("too small")
            continue

        # create a dataframe with winning vector at certain points in the section
        df = winners(vecs, x_vals, m0, m1, y0, dx, dx_y)

        # plot poincare section and save
        try:
            plot(df, vecs, c[i][j], j, n_squares, index, test=False)
        except Exception as error:
            print(error)
            continue
        df.to_csv(os.path.join(
            "results", f"{n_squares} - {index}", "df - " + str(j)), index=False)

        # make section object that define winning vector and line equations for boundaries of subsections
        sec_list = sec_setup(df, dx_y)
        secs = sec_comp(sec_list, dx)
        with open(os.path.join("results", f"{n_squares} - {index}", "secs - " + str(j) + ".dill"), 'wb') as f:
            dill.dump(secs, f)

        times = [1]
        if improved:
            times = time_comp(secs)

        # plot the pdf for each cusp
        pdf(list(df["time"]), times, dx*2, n_squares, index, j)

        # get covolume calculations
        covolume_list.append(covolume(secs))
        print("section done")

        break

with open(os.path.join("results", f"{n_squares} - {index}", "covolume.txt"), "w") as file:
    total = 0
    for item in covolume_list:
        total += item
        file.write(str(item) + "\n")
    file.write("\n" + str(total))
