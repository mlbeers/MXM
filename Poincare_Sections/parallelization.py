import unittest
from surface_dynamics.all import OrigamiDatabase, Origami
from sage.all import SymmetricGroup
from flatsurf import translation_surfaces
import numpy as np
from .Poincare import poincare_details, poincare_details, compute_winning_vecs_on_edges, setup, compute_winners_on_entire_section
from .utils import load_arrays_from_file  # testing
import time  # testing
from multiprocessing import Pool
import os  # testing
import concurrent.futures
from .computations import unit_horizontal_saddle, is_unobstructed, generate_permutations, generate_vectors, generate_poincare_section_details, compute_winners


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
    try:
        alphas, c_matrices, _, _, _, generators, eigenvectors = poincare_details(
            perm, vecs)

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ComputationsTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_parallel_generating_vectors_small(self):
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

    def test_parallel_generating_vectors_medium(self):
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

    def test_parallel_poincare_details_large(self):
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

    def test_parallel_poincare_details_small(self):
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
