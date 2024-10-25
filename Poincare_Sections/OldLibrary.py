import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import pwlf
import sympy as sym
from sympy import Symbol
from sympy import solve, lambdify
import os
import math
from surface_dynamics.all import *
from time import time
from scipy import integrate
import unittest
from surface_dynamics.all import Origami
from Library import DetailsError
from utils import load_arrays_from_file  # testing
from sage.all import matrix  # testing

# For each cusp of the square tiled surface, computes a generating matrix and
# its eigenvectors, then finds the shortest saddle connection in the direction # of the eigenvector by producing a big list of saddle connections and then checking all of them.
#
# inputs:
#   perm: a permutation that defines a square tiled surface
#   vecs0: list of saddle connections on this STS to check (this is large)
#
# outputs:
#   alphas:
#   Cs:
#   C_invs:
#   eigs: eigenvalues
#   Ms: C*generator*C^inv
#   generators:
#   eigenvecs: eigenvectors
##  old implementation of the Library.py/poincare_details

def poincare_details_old(perm, vecs0, generators):
    # find the eigenvectors for each generator and make sure they are 1
    eigs = []
    for matrix in generators:
        eig1, eig2 = matrix.eigenvalues()
        if eig1 == eig2:
            if eig1 == 1:
                eigs.append(eig1)
            else:
                raise ValueError("Eigenvalue not equal to 1")
        else:
            raise ValueError("Different eigenvalues")
    # find the eigenvectors for each generator
    eigenvecs = []
    for matrix in generators:
        vec = matrix.eigenvectors_right()[0][1][0]
        vec = np.array([[vec[0]], [vec[1]]])
        eigenvecs.append(vec)
    # find the magnitude, slope, x-direction, and y-direction of each eigenvector
    saddle_vecs = []
    for vec in eigenvecs:
        mag_vec = ((vec[0]**2 + vec[1]**2)**0.5)[0]
        if vec[0] == 0:
            slope_vec = float("inf")
        else:
            slope_vec = (vec[1]/vec[0])[0]

        if vec[0] >= 0:
            x_sign_vec = 1
        else:
            x_sign_vec = -1
        if vec[1] >= 0:
            y_sign_vec = 1
        else:
            y_sign_vec = -1

        saddle_vec = None
        check = 0

        # find the magnitude, slope, x-direction, and y-direction of each saddle connection
        for saddle in vecs0:
            mag_saddle = ((saddle[0]**2 + saddle[1]**2)**0.5)[0]
            if saddle[0] == 0:
                slope_saddle = float("inf")
            else:
                slope_saddle = (saddle[1]/saddle[0])[0]

            if saddle[0] >= 0:
                x_sign_saddle = 1
            else:
                x_sign_saddle = -1
            if saddle[1] >= 0:
                y_sign_saddle = 1
            else:
                y_sign_saddle = -1

            # find the smallest saddle connection that is in the same direction and has the same slope as the given eigenvector and add it to a list
            if slope_vec == slope_saddle:
                if x_sign_vec == x_sign_saddle:
                    if y_sign_vec == y_sign_saddle:
                        if check == 0:
                            saddle_vec = saddle
                            mag = mag_saddle
                            check += 1
                        elif mag_saddle < mag:
                            saddle_vec = saddle
                            mag = mag_saddle
        if check == 0:
            raise DetailsError(f"No saddle vec for eigenvector {vec}")
        saddle_vecs.append(saddle_vec)

    # find the counter-clockwise angle from the x-axis to the eigenvectors
    thetas = []
    for vec in saddle_vecs:
        mag = np.linalg.norm(vec)  # Calculate magnitude
        dot_product = np.dot([1, 0], vec)  # Dot product with (1,0)
        theta = np.arccos(dot_product / mag)  # Calculate angle
        if vec[1] < 0:
            theta = 2 * np.pi - theta
        thetas.append(theta)
    
    # find the rotation matrix that takes the vector (1,0) to the vector in the direction of each eigenvector
    rots = []
    for theta in thetas:
        rot = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]])
        rots.append(rot.reshape(2,2))
    
    # find c_inv and c
    Cs = []
    C_invs = []
    for i in range(len(rots)):
        mag = np.linalg.norm(saddle_vecs[i])
        if mag == 0:
            raise ValueError("Magnitude of saddle vector is zero.")
        S = np.array([[mag, 0], [0, 1/mag]])
        c_inv = rots[i] @ S
        c = np.linalg.inv(c_inv)
        Cs.append(c)
        C_invs.append(c_inv)

    epsilon = 0.005

    # alpha is the top right value of the matrix M = c @ generator @ c_inv. 
    # M must have 1s on the diagonal and 0 in the bottom left
    alphas = []
    Ms = []
    for i in range(len(generators)):
        M = Cs[i] @ generators[i] @ C_invs[i]
        Ms.append(M)
    
        # Check if the conditions are met
        if not (abs(M[0][0] - 1) < epsilon and 
                abs(M[1][1] - 1) < epsilon and 
                abs(M[1][0]) < epsilon):
            raise ValueError(
                f"Wrong conjugate matrix\nC: {Cs[i]}\nC_inv: {C_invs[i]}\nM: {M}\ngenerator: {generators[i]}\nsaddle: {saddle_vecs[i]}, \neigenvecs: {eigenvecs[i]}, \ndeterminant C: {np.linalg.det(Cs[i])}, \n {Cs[i]@saddle_vecs[i]}, \n {C_invs[i]@np.array([[1],[0]])}")
        
        alphas.append(round(M[0][1], 5))
        
    return alphas, Cs, C_invs, eigs, Ms, generators, eigenvecs
