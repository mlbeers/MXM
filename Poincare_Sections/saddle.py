
# sage saddle.py horizontal vertical 256R squares index


import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from time import time
#!/bin/env -S sage -python
from flatsurf import *
from sage.all import *
import sys
import pandas as pd
import re
import os


def vectors(h, v, length=256):
    from flatsurf import translation_surfaces
    S = SymmetricGroup(len(h))
    T = translation_surfaces.origami(S(h), S(v))
    T = T.erase_marked_points()
    from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
    TT = to_pyflatsurf(T)
    from pyflatsurf import flatsurf
    C = TT.connections().bound(flatsurf.Bound(length))
    connections = [
        vector((ZZ(str(c.vector().x())), ZZ(str(c.vector().y())))) for c in C]
    vecs = []
    for vec in connections:
        item = np.array([[vec[0]], [vec[1]]])
        vecs.append(item)
    print(len(vecs))
    return vecs


h = str(sys.argv[1])
v = str(sys.argv[2])
n = int(sys.argv[3])
squares = int(sys.argv[4])
index = int(sys.argv[5])


name = "vecs" + str(squares) + "-" + str(index) + ".npy"
vecs0 = vectors(h, v, n)


def save_arrays_to_file(file_path, arrays_list):
    # Save arrays to a single NumPy file
    np.save(file_path, arrays_list)


save_arrays_to_file(os.path.join("vecs", name), vecs0)
