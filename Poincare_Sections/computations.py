import unittest
from surface_dynamics.all import OrigamiDatabase, Origami
from sage.all import SymmetricGroup
from flatsurf import translation_surfaces
import numpy as np
from .Poincare import poincare_details, poincare_details, winners
from .utils import load_arrays_from_file, save_arrays_to_file  # testing
import time  # testing
from multiprocessing import Pool
import os  # testing
import concurrent.futures


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
        # direction = sc.direction  # TODO: what's this for?
        if vec not in vectors:
            if (vec[0] >= -length/20 and vec[0] <= length/20) and (vec[1] >= -length/20 and vec[1] <= length/20):
                vectors.append(vec)

    vectors = [np.array([[v[0]], [v[1]]]) for v in vectors]
    return vectors


# mildly faster. the bottleneck is in generating saddle connections
# TODO: would be even slightly faster without all the casting at the end.

def parallel_generate_vectors(permutation, radius=200):
    h, v = str(permutation).split("\n")
    S = SymmetricGroup(len(h))
    T = translation_surfaces.origami(S(h), S(v))
    # T = T.erase_marked_points() # only works if pyflatsurf is installed
    saddle_connections = T.saddle_connections(
        radius)  # TODO: is this only first quadrant?

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        vectors = executor.map(inner_generate_vectors, [
            (sc, radius) for sc in saddle_connections])

    # remove None (when inner_generate_vectors doesn't find a vec within radius)
    vectors = [v for v in vectors if v is not None]
    # cast to hashable type for deduplicating
    vectors = [(v[0], v[1]) for v in vectors]
    # deduplicate, note: this reorders the list
    vectors = list(set(vectors))
    # cast to type needed for later computations
    vectors = [np.array([[v[0]], [v[1]]]) for v in vectors]

    return vectors


def inner_generate_vectors(sc_data):
    saddle_connection, radius = sc_data
    vec = saddle_connection.holonomy().n()
    # direction = saddle_connection.direction  # TODO: what's this?
    if (vec[0] >= -radius/20 and vec[0] <= radius/20) and (vec[1] >= -radius/20 and vec[1] <= radius/20):
        return vec


# wraps the poincare_details function in a try catch because
# when the number of vectors used in this simulation is too
# low, the chances of there not being a precomputed vector in
# the direction of the saddle connection is reasonably high.
# this function can take a while!
# - permutation: a permutation defining a STS
# - vectors: pregenerated list of vectors in the direction of saddle connections on this STS

def try_poincare_details(sts_data):
    perm, vecs = sts_data
    details = []
    for i in range(1):
        try:
            alphas, c_matrices, _, _, _, generators, eigenvectors = poincare_details(
                (perm, vecs))

            details.append((alphas, c_matrices, generators, eigenvectors))
        except:
            pass

    return details


def init_pool_processes(p, v):
    global permutation
    global vectors
    permutation = p
    vectors = v


def threaded_generate_poincare_section_details(sts_data, trys):
    details = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        details = executor.map(try_poincare_details, [
            sts_data for i in range(trys)])

    details = [x for x in details if not len(x) == 0]
    return details


def multiprocess_generate_poincare_section_details(sts_data, trys):
    details = []
    with Pool() as p:
        details = p.map(try_poincare_details, [
            sts_data for i in range(trys)])

    details = [x for x in details if not len(x) == 0]
    return details


def generate_poincare_section_details(sts_data, trys):
    permutation, vectors = sts_data

    details = []
    for i in range(trys):
        try:
            alphas, c_matrices, _, _, _, generators, eigenvectors = poincare_details(
                (permutation, vectors))

            details.append((alphas, c_matrices, generators, eigenvectors))
        except:
            pass

    return details


# Run computations for individual cusps
# j represents the jth cusps
# change i to get new output for the jth cusp
# output:
# - number of vectors
# - the x and y coords of the "section vector"
# - poincare section plot

def compute_on_cusp(details, vectors):
    (alphas, c_matrices, eigenvectors, _) = details
    # 0, 24, 24, 0, 6, 13, 0, 3, 3, 1
    i = 24  # depends on number of times you run poincare_details
    j = 1  # depends on number of cusps of STS
    n_squares = 7
    dx = 0.0005
    index = 4
    vecs, x_vals, m0, m1, x0, y0, dx_y = setup(
        alphas[i][j], c_matrices[i][j], eigenvectors[i][j], vectors, dx, False)
    df = winners(vecs, x_vals, m0, m1, y0, dx, dx_y)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    def test_generating_vectors_small(self):
        permutations = generate_permutations(7)
        start_1 = time.time()
        vectors = generate_vectors(permutations[4], 200)
        end_1 = time.time()
        self.assertEqual(len(vectors), 256)  # TODO: why is this 256?

        start_2 = time.time()
        vectors = parallel_generate_vectors(permutations[4], 200)
        end_2 = time.time()
        self.assertEqual(len(vectors), 256)  # TODO: why is this 256?
        print(
            f'parallel speed: {end_2-start_2}\n regular speed: {end_1-start_1}')
        self.assertGreater(end_2-start_2, end_1-start_1)

    def test_generating_vectors_medium(self):
        permutations = generate_permutations(7)
        start_1 = time.time()
        vectors = generate_vectors(permutations[4], 1000)
        end_1 = time.time()
        # self.assertEqual(len(vectors), 256)  # TODO: why is this 256?

        start_2 = time.time()
        vectors = parallel_generate_vectors(permutations[4], 1000)
        end_2 = time.time()
        # self.assertEqual(len(vectors), 256)  # TODO: why is this 256?
        print(
            f'parallel speed: {end_2-start_2}\n regular speed: {end_1-start_1}')
        self.assertGreater(end_2-start_2, end_1-start_1)

    def test_poincare_details_large_data(self):
        # mock data
        # this is permutations[3] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5)(6,7)', '(1,2,3,4,6)(5,7)')
        # this is corresponding vector data computed with full library
        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-3.npy"))
        self.assertEqual(len(vecs), 907654)

        # print('running threaded poincare section computation')
        # t0 = time.time()
        # details_2 = threaded_generate_poincare_section_details(
        #     (perm, vecs), 10)
        # t1 = time.time()
        # print(f'time: {t1-t0}')
        # print(details_2)

        print('running multiprocess poincare section computation')
        t0 = time.time()
        details_3 = multiprocess_generate_poincare_section_details(
            (perm, vecs), 3)
        t1 = time.time()
        print(f'time: {t1-t0}')
        print(details_3)

        print('running regular poincare section computation')
        t0 = time.time()
        details_1 = generate_poincare_section_details((perm, vecs), 3)
        t1 = time.time()
        print(f'time: {t1-t0}')
        print(details_1)

    def test_poincare_details_small_data(self):
        # this is permutations[4] from generate_permutations(7)
        print('generating test data')
        perm = Origami(
            '(1)(2)(3)(4,5,6,7)', '(1,2,3,4)(5)(6)(7)')

        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "test_data_4.npy"))
        self.assertEqual(len(vecs), 1912)

        print('running threaded poincare section computation')
        t0 = time.time()
        details_2 = threaded_generate_poincare_section_details((perm, vecs), 1)
        t1 = time.time()
        print(f'time: {t1-t0}')
        print(details_2)

        print('running multiprocess poincare section computation')
        t0 = time.time()
        details_3 = multiprocess_generate_poincare_section_details(
            (perm, vecs), 1)
        t1 = time.time()
        print(f'time: {t1-t0}')
        print(details_3)

        print('running regular poincare section computation')
        t0 = time.time()
        details_1 = generate_poincare_section_details((perm, vecs), 1)
        t1 = time.time()
        print(f'time: {t1-t0}')
        print(details_1)
        # output = compute_on_cusp(details, vecs)

    def test_winners(self):
        # this is permutations[3] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5)(6,7)', '(1,2,3,4,6)(5,7)')
        # this is corresponding vector data computed with full library
        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-3.npy"))
        details = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "poincare_details_output.txt"))

        compute_on_cusp(details, vecs)
