import numpy as np
from matplotlib import pyplot as plt
from Library import *
from Library import Section
from vector import generate_vectors
from time import time
from flatsurf import *
from sage.all import *
import sys
import os

t0 = time()

# number of squares for STS
n_squares = int(sys.argv[1])
# index to start at
index = int(sys.argv[2])

# list of permuatations
permutations = perms_list(n_squares)
# generate saddle connections on surface and save file
name = "vecs" + str(n_squares) + "-" + str(index) + ".npy"

vecs0 = generate_vectors(permutations[index], 2000)
save_arrays_to_file(os.path.join("vecs", name), vecs0)

t1 = time()
print(f"vectors done: {(t1-t0)/60**2}\n")
print("-----------------------------------------------------------------------------------")