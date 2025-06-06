import numpy as np
#from matplotlib import pyplot as plt
from Library import *
#from Library import Section
from time import time
from flatsurf import *
from sage.all import *
#import sys
#import pandas as pd
#import re
#import os
import time  # testing


# generate vectors for saddle connections on Square Tiled
# Surfaces (STS)
# - perm: a permutation defining a STS
# - length: TODO: what is this parameter?
def generate_vectors(perm, length=256):
    a = str(perm)
    h, v = a.split("\n")
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

# generate vectors for saddle connections on STS, alternate method
def vectors2(perm, length=200):
    a = str(perm)
    h, v = a.split("\n")
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
        item = np.array([[vec[0]], [vec[1]]])
        vecs.append(item)
    return vecs


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ComputationsTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_generating_vectors_small(self):
        print(f'\nTesting: generate_vectors() with small data')
        perm = '(1)(2)(3)(4,5,6,7)\n(1,2,3,4,7,6,5)'
        t0 = time.time()
        vectors = generate_vectors(perm, 200)
        t1 = time.time()
        self.assertEqual(len(vectors), 256)  # TODO: why is this 256?
        print(f'  runtime: {t1-t0}')

    def test_generating_vectors_medium(self):
        print(f'\nTesting: generate_vectors() with medium data')
        perm = Origami(
            '(1)(2)(3)(4,5,6,7)', '(1,2,3,4,7,6,5)')
        t0 = time.time()
        vectors = generate_vectors(perm, 1000)
        t1 = time.time()
        # self.assertEqual(len(vectors), 256)  # TODO: why is this 256?
        print(f'  runtime: {t1-t0}')
