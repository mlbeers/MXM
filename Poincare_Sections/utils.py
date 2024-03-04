import os
import numpy as np
import unittest  # for testing
import logging  # for testing

# take files from saddle.py and load them into notebook


def load_arrays_from_file(file_path):
    # Load arrays from the NumPy file
    arrays_list = np.load(file_path, allow_pickle=True)

    # Ensure each element in the list is a NumPy array
    arrays_list = [np.array(array) for array in arrays_list]

    return arrays_list


class UtilsTestSuite(unittest.TestCase):

    """Test loading vectors with two different ways of constructing the correct relative paths"""

    def test_generating_permutations(self):
        # print(os.path.abspath(".")) # use to check base path
        vecs = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-3.npy"))
        self.assertEqual(len(vecs), 907654)
        vecs = load_arrays_from_file("Poincare_sections/vecs/vecs7-4.npy")
        self.assertEqual(len(vecs), 3618904)
