#sage saddle.py horizontal vertical 500000 squares index

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

def vectors(h, v, length = 200):
    S = SymmetricGroup(len(h))
    T = translation_surfaces.origami(S(h), S(v))
    T = T.erase_marked_points()
    sc_list = T.saddle_connections(length)
    slopes_all = []
    for item in sc_list:
        vec = item.holonomy().n()
        direction = item.direction
        if vec not in slopes_all:
            if vec[0] >= -length/20 and vec[0] <= length/20:
                if vec[1] >= -length/20 and vec[1] <= length/20:
                    slopes_all.append(item.holonomy().n())         
    vecs = []
    for vec in slopes_all:
        item = np.array([[vec[0]],[vec[1]]])
        vecs.append(item)
    return vecs


h = str(sys.argv[1])
v = str(sys.argv[2])
n = int(sys.argv[3])
squares = int(sys.argv[4])
index = int(sys.argv[5])


name = "vecs" + str(squares) + "-" + str(index) + ".txt"
vecs0 = vectors(h,v,n)

f = open(os.path.join("vecs", name),'w')
f.write(str(vecs0))
f.close()