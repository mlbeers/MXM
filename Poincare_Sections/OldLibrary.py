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

# old implementation of the Library.py/poincare_details
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

# old implementation of the Library.py/winners
# It first computes winning vectors from all possible vectors on points on the
# bottom and vertical edge of the poincare section. Then it computes winners
# for all points on the poincare section from the set of winners from the sides.
#
# inputs:
#   vecs: array of vectors, after being multiplied by C matrix
#   x_vals: array of x-values in the poincare section, all between [0,1]
#   m0: slope of the top line of the poincare section
#   m1: slope of the bottom line of the poincare section
#   y0: 1/y0 is the y-coordinate of the top left corner point of the section
#   dx: step size for separation of points to sample in the section
#   dx_y: step size for separation of points to sample in the section
#
# output: TODO
def winners_old(vecs, x_vals, m0, m1, y0, dx, dx_y):
    # dictionary for plotting
    saddle_dict = {}
    saddle_dict["x"] = []
    saddle_dict["y"] = []
    saddle_dict["lab"] = []
    saddle_dict["vec"] = []
    saddle_dict["time"] = []
    possible_vecs = []
    winners = []

    t0 = time()
    # top edge
    dz = 1/100
    for a in np.arange(dz, 1, dz):
        check = 0
        winner_slope = None
        winner = None
        Mab = np.array([[a, 1-dx/y0], [0, a]])
        for vec in vecs:
            new = Mab@vec
            if float(new[0][0]) == 0:
                continue
            x = float(new[0][0])
            y = float(new[1][0])
            if y/x <= 0:
                continue
            if x <= 1 and x > 0:
                if winner_slope == None:
                    winner_slope = y/x
                    winner = vec
                    continue
        # if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                elif abs(y/x - winner_slope) <= dx/1000:
                    if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                        winner = vec
                        continue
                elif y/x < winner_slope:
                    winner_slope = y/x
                    winner = vec
                    continue
        winners.append(winner)
    t1 = time()
    print("top done: " + str(t1 - t0))

    # could re-include if we need to do computations on all sides
    # # diagonal
    # for a in np.arange(0 + dz, 1, dz):
    #     check = 0
    #     winner_slope = None
    #     winner = None
    #     Mab = np.array([[a, m1*a + 1-dx/y0 + dx_y], [0, a]])
    #     for vec in vecs:
    #         new = Mab@vec
    #         if float(new[0][0]) == 0:
    #             continue
    #         x = float(new[0][0])
    #         y = float(new[1][0])
    #         if y/x <= 0:
    #             continue
    #         if x <= 1 and x > 0:
    #             if winner_slope == None:
    #                 winner_slope = y/x
    #                 winner = vec
    #                 continue
    #    # if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
    #             elif abs(y/x - winner_slope) <= dx/1000:
    #                 if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
    #                     winner = vec
    #                     continue
    #             elif y/x < winner_slope:
    #                 winner_slope = y/x
    #                 winner = vec
    #                 continue
    #     winners.append(winner)

    # side edge
    t0 = time()
    y_vals = np.arange(m1 + (1-dx)/y0 + dx_y, m0 +
                       (1-dx)/y0 - dx_y, dz*(m0-m1))
    for b in y_vals:
        check = 0
        winner_slope = None
        winner = None
        Mab = np.array([[1 - dx, b], [0, 1-dx]])
        for vec in vecs:
            new = Mab@vec
            if float(new[0][0]) == 0:
                continue
            x = float(new[0][0])
            y = float(new[1][0])
            if y/x <= 0:
                continue
            if x <= 1 and x > 0:
                if winner_slope == None:
                    winner_slope = y/x
                    winner = vec
                    continue
        # if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                elif abs(y/x - winner_slope) <= dx/1000:
                    if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                        winner = vec
                        continue
                elif y/x < winner_slope:
                    winner_slope = y/x
                    winner = vec
                    continue
        winners.append(winner)
    t1 = time()
    print("side done: " + str(t1 - t0))

    winners2 = []
    for winner in winners:
        try:
            if winner == None:
                continue
        except:
            winners2.append(winner)
    possible_vecs1 = np.unique(winners2, axis=0)
    for item in possible_vecs1:
        possible_vecs.append(item)

    global label_dict
    label_dict = {}

    # dictionary for vector labels
    for i in range(len(possible_vecs)):
        label_dict[i] = possible_vecs[i]

    # for each vector, there is a time function defined as f(a,b) where a,b are points in the poincare section
    global t_dict
    t_dict = {}

    x, y, t = sym.symbols('x y t')
    Mab = np.array([[x, y], [0, 1/x]])
    horo = np.array([[1, 0], [-t, 1]])

    for i in range(len(possible_vecs)):
        # apply Mab matrix, perform horocycle flow and find time t to horizontal
        a = horo@(Mab@possible_vecs[i])
        t_dict[i] = lambdify([x, y], solve(a[1][0], t)[0])

    # for each point (a,b) in the poincare section, apply the Mab matrix to each vector and look for "winners". Winners have smallest possible slope that is greater than zero and 0 < x-component <= 1
    for a in x_vals:
        y_vals = np.arange(m1*a + 1/y0 + dx_y, m0*a + 1/y0 - dx_y, dx_y)
        for b in y_vals:
            check = 0
            winner_slope = None
            winner = None
            Mab = np.array([[a, b], [0, 1/a]])
            for vec in possible_vecs:
                new = Mab@vec
                if float(new[0][0]) == 0:
                    continue
                x = float(new[0][0])
                y = float(new[1][0])
                if y/x <= 0:
                    continue
                if x <= 1 and x > 0:
                    if winner_slope == None:
                        winner_slope = y/x
                        winner = vec
                        continue
            # if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                    elif abs(y/x - winner_slope) <= dx/1000:
                        if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                            winner = vec
                            continue
                    elif y/x < winner_slope:
                        winner_slope = y/x
                        winner = vec
                        continue

            saddle_dict["x"].append(a)
            saddle_dict["y"].append(b)
            saddle_dict["vec"].append(winner)
            for i in range(len(possible_vecs)):
                if np.array_equal(winner, possible_vecs[i]):
                    check += 1
                    saddle_dict["lab"].append(i)
                    saddle_dict["time"].append(t_dict[i](a, b))
                    if saddle_dict["time"][-1] < 0:
                        saddle_dict["time"][-1] = 1000
            # if there is no winner, at (a,b), add a label so the df and plot can still be made. These section will later be made blank for trouble-shoooting
            if check == 0:
                saddle_dict["lab"].append(len(vecs))
                saddle_dict["time"].append(1000)

    df = pd.DataFrame.from_dict(saddle_dict)
    return df
