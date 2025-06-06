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

vec_length = int(sys.argv[1])
folder = str(sys.argv[2])

with open(os.path.join("results", folder, "perm.dill"), 'rb') as f:
    perm = dill.load(f)
    
name = "vecs" + folder + ".npy"

vecs0 = generate_vectors(perm, vec_length)
save_arrays_to_file(os.path.join("vecs", name), vecs0)

t1 = time()
print(f"vectors done: {(t1-t0)/60**2}\n")
print("-----------------------------------------------------------------------------------")