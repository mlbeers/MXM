import numpy as np
from matplotlib import pyplot as plt
from Library import *
from Library import Section
from .vector import generate_vectors
import pandas as pd
from time import time
from flatsurf import *
from sage.all import *
import sys
import pandas as pd
import re
import os


# number of squares for STS
n_squares = int(sys.argv[1])
# index to start at
start_index = int(sys.argv[2])
# number of squares to skip count by if we are running the script in multiple tmux windows
increment = int(sys.argv[3])

# list of permuatations
permutations = perms_list(n_squares)

for i in range(start_index, len(permutations), increment):
    # get the horizontal and vertical permutations
    h = str(permutations[i].r())
    v = str(permutations[i].u())
    for num in range(1, n_squares+1):
        if str(num) not in h:
            h += "(" + str(num) + ")"
        if str(num) not in v:
            v += "(" + str(num) + ")"
    # generate saddle connections on surface and save file
    name = "vecs" + str(n_squares) + "-" + str(i) + ".npy"

    vecs0 = generate_vectors(h, v, 2000)
    save_arrays_to_file(os.path.join("vecs", name), vecs0)
