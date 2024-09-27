import unittest
from surface_dynamics.all import OrigamiDatabase, Origami
from sage.all import SymmetricGroup
from flatsurf import translation_surfaces
import numpy as np
# , compute_winning_vecs_on_edges, compute_winners_on_entire_section
from .Library import poincare_details, poincare_details, setup
from .utils import load_arrays_from_file  # testing
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

def compute_winners(details, vectors):
    (alphas, c_matrices, eigenvectors, _) = details
    # 0, 24, 24, 0, 6, 13, 0, 3, 3, 1
    i = 0  # depends on number of times you run poincare_details
    j = 1  # depends on number of cusps of STS
    n_squares = 7
    dx = 0.0005
    index = 4
    vecs, x_vals, m0, m1, x0, y0, dx_y = setup(
        alphas[i], c_matrices[i], eigenvectors[i], vectors, dx, False)

    edge_winners = compute_winning_vecs_on_edges(
        vecs, x_vals, m0, m1, y0, dx, dx_y)
    winners_data_frame = compute_winners_on_entire_section(
        edge_winners, x_vals, m0, m1, y0, dx, dx_y)
    return winners_data_frame


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ComputationsTestSuite(unittest.TestCase):
    """Basic test cases."""

    # TODO: why these numbers?
    def test_generating_permutations(self):
        print(f'Testing: generate_permutations()')
        permutations = generate_permutations(5)
        self.assertEqual(len(permutations), 0)
        permutations = generate_permutations(6)
        self.assertEqual(len(permutations), 37)
        permutations = generate_permutations(7)
        self.assertEqual(len(permutations), 92)
        permutations = generate_permutations(8)
        self.assertEqual(len(permutations), 352)

    # TODO: Not all permutuations work?

    def test_poincare_details_large_data(self):
        print(f'\nTesting: poincare_details() with large data')
        # Set up mock data
        # this is permutations[3] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5)(6,7)', '(1,2,3,4,6)(5,7)')
        # this is corresponding vector data computed with full library
        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-3.npy"))
        self.assertEqual(len(vecs), 907654)

        # Run test
        t0 = time.time()
        details_1 = generate_poincare_section_details((perm, vecs), 3)
        t1 = time.time()
        print(f'  runtime: {t1-t0}')
        # print(details_1)

    def test_poincare_details_small_data(self):
        print(f'\nTesting: poincare_details() with small data')
        # Set up mock data
        # this is permutations[4] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5,6,7)', '(1,2,3,4)(5)(6)(7)')

        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "test_data_4.npy"))
        self.assertEqual(len(vecs), 1912)

        # Run test
        t0 = time.time()
        details_1 = generate_poincare_section_details((perm, vecs), 1)
        t1 = time.time()
        print(f'  runtime: {t1-t0}')
        # print(details_1)
        # output = compute_on_cusp(details, vecs)