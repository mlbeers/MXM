import numpy as np
from matplotlib import pyplot as plt
from Library import *
from Library import Section
import pandas as pd
from time import time
from flatsurf import *
from sage.all import *
import sys
import pandas as pd
import re
import os

def vectors(h, v, length = 256):
    from flatsurf import translation_surfaces
    S = SymmetricGroup(len(h))
    T = translation_surfaces.origami(S(h), S(v))
    T = T.erase_marked_points()
    from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
    TT = to_pyflatsurf(T)
    from pyflatsurf import flatsurf
    C = TT.connections().bound(flatsurf.Bound(length))
    connections = [vector((ZZ(str(c.vector().x())), ZZ(str(c.vector().y())))) for c in C]
    vecs = []
    for vec in connections:
        item = np.array([[vec[0]],[vec[1]]])
        vecs.append(item)
    print(len(vecs))
    return vecs
    
#number of squares for STS
n_squares = int(sys.argv[1])
#index to start at
start_index = int(sys.argv[2])
#number of squares to skip count by if we are running the script in multiple tmux windows
increment = int(sys.argv[3])

#list of permuatations
permutations = perms_list(n_squares)

for i in range(start_index, len(permutations), increment):
    #get the horizontal and vertical permutations
    h = str(permutations[i].r())
    v = str(permutations[i].u())
    for num in range(1, n_squares+1):
        if str(num) not in h:
            h += "(" + str(num) + ")"
        if str(num) not in v:
            v += "(" + str(num) + ")"
    #generate saddle connections on surface and save file
    name = "vecs" + str(n_squares) + "-" + str(i) + ".npy"
    vecs0 = vectors(h, v, 2000)
    save_arrays_to_file(os.path.join("vecs", name), vecs0)