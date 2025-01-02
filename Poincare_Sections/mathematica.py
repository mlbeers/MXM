import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import os
import pwlf
import os
from surface_dynamics.all import *
from Library import *
from Library import Section
import math
from time import time
import copy
from scipy import integrate
import sympy as sym
from sympy import Symbol
from sympy import solve, lambdify
import traceback
import dill
import sys
import unittest
from surface_dynamics.all import Origami
from utils import load_arrays_from_file  # testing

from sage.all import *
import numpy as np
from fractions import Fraction as frac
M = mathematica

# num = 1
# matrix_string = f"A = {{{{{-344 * num}, {225 * num}}}, {{{-529 * num}, {346 * num}}}}}"
#print(matrix_string)

# A = M(matrix_string)  # Pass the constructed string to Mathematica
# eigenvalues = M('eigenvalues = Eigenvalues[A]')
# eig1 = int(M('eigenvalues[[1]]'))
# eigenvectors = M('eigenvectors = Eigenvectors[A]')
# scaledEigenvector = M('scaledEigenvectors = Map[LCM @@ Denominator[#]*# &, eigenvectors][[1]]')
# S0, J = M('{S0, J} = JordanDecomposition[A]')
# S = M('S = Inverse[S0]')

def poincare_details(perm, vecs0, generators):
    M = mathematica
    # find the eigenvectors for each generator and make sure they are 1
    eigs = []
    eigenvecs = []
    for matrix in generators:
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]

        matrix_string = f"A = {{{{{a}, {b}}}, {{{c}, {d}}}}}"
        A = M(matrix_string)
        eigenvalues = M('eigenvalues = Eigenvalues[A]')
        eig1 = int(M('eigenvalues[[1]]'))
        eig2 = int(M('eigenvalues[[2]]'))

        if eig1 == eig2:
            if eig1 == 1:
                eigs.append(eig1)
            else:
                raise ValueError("Eigenvalue not equal to 1")
        else:
            raise ValueError("Different eigenvalues")
        
        eigenvectors = M('eigenvectors = Eigenvectors[A]')
        scaledEigenvector = M('scaledEigenvectors = Map[LCM @@ Denominator[#]*# &, eigenvectors][[1]]')
        scaled_x = int(scaledEigenvector.sage()[0])
        scaled_y = int(scaledEigenvector.sage()[1])
        vec = np.array([[scaled_x], [scaled_y]])
        eigenvecs.append(vec)

    print(eigenvecs)
    # find the magnitude, slope, x-direction, and y-direction of each eigenvector
    def get_magnitude_slope_sign(vectors):
        # Stack the list of 2D numpy arrays into a 2D array
        vectors = np.hstack(vectors)  # Concatenate into a 2D array
        
        magnitudes = np.linalg.norm(vectors, axis=0)
        
        # Make sure the output array for division is explicitly float
        slopes = np.divide(vectors[1], vectors[0], out=np.full_like(vectors[1], float('inf'), dtype=np.float64), where=vectors[0] != 0)
        
        x_signs = np.sign(vectors[0])
        y_signs = np.sign(vectors[1])
        
        return magnitudes, slopes, x_signs, y_signs

    eigen_mags, eigen_slopes, eigen_x_signs, eigen_y_signs = get_magnitude_slope_sign(eigenvecs)
    saddle_mags, saddle_slopes, saddle_x_signs, saddle_y_signs = get_magnitude_slope_sign(vecs0)

    # Find corresponding saddle vectors for eigenvectors
    saddle_vecs = []
    for i in range(len(eigenvecs)):
        slope_vec = eigen_slopes[i]
        x_sign_vec = eigen_x_signs[i]
        y_sign_vec = eigen_y_signs[i]

        # Identify matches by slope, x, and y direction
        slope_matches = (slope_vec == saddle_slopes)
        x_sign_matches = x_sign_vec == saddle_x_signs
        y_sign_matches = y_sign_vec == saddle_y_signs

        valid_saddles = np.where(slope_matches & x_sign_matches & y_sign_matches)[0]

        if len(valid_saddles) == 0:
            raise DetailsError(f"No saddle vec for eigenvector {eigenvecs[i]}")

        # Select saddle with smallest magnitude
        smallest_idx = valid_saddles[np.argmin(saddle_mags[valid_saddles])]
        saddle_vecs.append(vecs0[smallest_idx])

    saddle_vecs = np.array(saddle_vecs)
    print(saddle_vecs)
    
    # find c_inv and c
    Cs = []
    alphas = []
    Ms = []
    Ss = []
    for matrix, saddle_vec, eigen_vec in zip(generators, saddle_vecs, eigenvecs):
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]
    
        matrix_string = f"A = {{{{{a}, {b}}}, {{{c}, {d}}}}}"
        A = M(matrix_string)
    
        S0, J = M('{S0, J} = JordanDecomposition[A]')
        S = M('S = Inverse[S0]')
    
        S_np = np.array(S.sage())
        a_ = frac(str(S_np[0][0]))
        b_ = frac(str(S_np[0][1]))
        c_ = frac(str(S_np[1][0]))
        d_ = frac(str(S_np[1][1]))
        c0 = np.array([[a_, b_], [c_, d_]], dtype=object)
        print(c0)
    
        M_np = np.array(J.sage())
        a_1 = int(M_np[0][0])
        b_1 = int(M_np[0][1])
        c_1 = int(M_np[1][0])
        d_1 = int(M_np[1][1])
        J = np.array([[a_1, b_1], [c_1, d_1]])
        print(J)
    
    
        if a_1 != 1 or c_1 != 0 or d_1 != 1:
            print()
            raise ValueError(f"wrong M: {J}")
            
        if saddle_vec[0][0] != 0:
            entry = saddle_vec[0][0] / eigen_vec[0][0]
        elif saddle_vec[1][0] != 0:
            entry = saddle_vec[1][0] / eigen_vec[1][0]
        else:
            raise ValueError("Both coordinates in saddle_vec are zero. Division by zero is not defined.")
    
        S = np.array([[entry, 0], [0, frac(str(1/entry))]], dtype=object)
        c = c0 @ S
        Cs.append(c)
        Ms.append(J)
        Ss.append(S)
        
        alphas.append(J[0][1])
            
        return alphas, Cs, Ss, eigs, Ms, generators, eigenvecs

################# TEST #####################
n_squares = 7
index = 0
permutations = perms_list(n_squares)
perm = permutations[index]

vec_file = "vecs" + str(n_squares) + "-" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))

gs = generators(perm, vecs0)
alphas, Cs, Ss, eigs, Ms, generators, eigenvecs = poincare_details(perm, vecs0, gs)

print(alphas)
print(eigs)
print(eigenvecs)
print(Cs)

