import unittest
from surface_dynamics.all import OrigamiDatabase, Origami
from sage.all import SymmetricGroup
from flatsurf import translation_surfaces
import numpy as np
from .Poincare import poincare_details
from .utils import load_arrays_from_file  # testing
import time  # testing
from multiprocessing import Pool
import os  # testing


# detect if a given orgiami has TODO: a unit horizontal
# saddle connection?

def unit_horizontal_saddle(origami):
    count = 0
    for vert in origami.vertices():
        tup = vert.up_right_tuple()
        for i in tup:
            for vert2 in origami.vertices():
                tup2 = vert2.up_right_tuple()
                if origami.r()(i) in tup2:
                    return True
    return False


# detect if a given origami is TODO: unobstructed ?

def is_unobstructed(origami):
    cusp_reps = origami.teichmueller_curve().cusp_representatives()

    for representative in cusp_reps:
        if not unit_horizontal_saddle(representative[0]):
            return False
    return True


# generate a list of permutations that define a Square
# Tiled Surface TODO: w/ some mysterious unobstructedness
# conditions ?
# - n: number of squares to use in generating the STSs
# output:
# - list of permutations: strings formatted like (..)\n(..)

def generate_permutations(n, **kwargs):
    obstructed = []
    Database = OrigamiDatabase()
    query_results = Database.query(nb_squares=n, **kwargs)

    for row in query_results:
        if not is_unobstructed(row):
            obstructed.append(row)
            for permutation in row.teichmueller_curve():
                obstructed.append(permutation)
    return obstructed


# generate vectors for saddle connections on Square Tiled
# Surfaces (STS)
# - perm: a permutation defining a STS
# - length: TODO: what is this parameter?

def generate_vectors(permutation, length=200):
    h, v = str(permutation).split("\n")
    S = SymmetricGroup(len(h))
    T = translation_surfaces.origami(S(h), S(v))
    # T = T.erase_marked_points() # only works if pyflatsurf is installed
    saddle_connections = T.saddle_connections(length)
    vectors = []
    for sc in saddle_connections:
        vec = sc.holonomy().n()
        direction = sc.direction  # TODO: what's this for?
        if vec not in vectors:
            if (vec[0] >= -length/20 and vec[0] <= length/20) and (vec[1] >= -length/20 and vec[1] <= length/20):
                vectors.append(vec)

    return [np.array([[v[0]], [v[1]]]) for v in vectors]

# wraps the poincare_details function in a try catch because
# when the number of vectors used in this simulation is too
# low, the chances of there not being a precomputed vector in
# the direction of the saddle connection is reasonably high.
# this function can take a while!
# - permutation: a permutation defining a STS
# - vectors: pregenerated list of vectors in the direction of saddle connections on this STS


def try_poincare_details(sts_data):
    permutation, vectors = sts_data
    details = []
    try:
        alphas, c_matrices, _, _, _, generators, eigenvectors = poincare_details(
            permutation, vectors)
        details.append((alphas, c_matrices, generators, eigenvectors))
    except:
        pass
    return details


def poincare_details_wrapper(sts_data, trys):
    details = []
    with Pool(20) as p:
        print(p.map(try_poincare_details, [sts_data for i in range(trys)]))


# ~~ Tests ~~

class ComputationsTestSuite(unittest.TestCase):
    """Basic test cases."""

    # TODO: why these numbers?
    def test_generating_permutations(self):
        permutations = generate_permutations(5)
        self.assertEqual(len(permutations), 0)
        permutations = generate_permutations(6)
        self.assertEqual(len(permutations), 37)
        permutations = generate_permutations(7)
        self.assertEqual(len(permutations), 92)
        permutations = generate_permutations(8)
        self.assertEqual(len(permutations), 352)

    # TODO: Not all permutuations work?

    def test_generating_vectors(self):
        permutations = generate_permutations(7)
        vectors = generate_vectors(permutations[4], 200)
        self.assertEqual(len(vectors), 200)  # TODO: why is this 256?

    def test_poincare_details_large_data(self):
        start = time.time()
        # mock test data
        # this is permutations[4] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5,6,7)', '(1,2,3,4)(5)(6)(7)')
        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-4.npy"))
        print(len(vecs))
        output = poincare_details_wrapper((perm, vecs), 600)
        end = time.time()
        print(end-start)
        print(output)

    def test_poincare_details_small_data(self):
        start = time.time()
        # this is permutations[4] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5,6,7)', '(1,2,3,4)(5)(6)(7)')
        vecs = generate_vectors(perm, 500)
        output = poincare_details_wrapper((perm, vecs), 600)
        end = time.time()
        print(end-start)
        print(output)
