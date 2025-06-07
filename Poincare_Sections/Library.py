from sage.all import *
from sage.all import matrix  # testing
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import pwlf
import os
from surface_dynamics.all import *
from surface_dynamics.all import Origami
from time import time
from scipy import integrate
import unittest
from utils import load_arrays_from_file  # testing
from fractions import Fraction as frac
import sympy as sym
from sympy import Matrix, sqrt, solve, lambdify, Eq, Symbol
import dill
from collections import defaultdict
from estimated_library import *

def compute_poincare_sections(vecs0, a, c, e, dx, dy, dx_frac, dy_frac, j, folder, estimated_output=False):
    print(f"id: {os.getpid()}")
    list_as = []
    list_cs = []
    list_es = []
    for a_, c_, e_ in zip(a, c, e):
        list_as.append(a_[j])
        list_cs.append(c_[j])
        list_es.append(e_[j])

    # Sort a, c, and e together based on the absolute value of the lower-right + lower-left entry of c
    sorted_pairs = sorted(zip(list_as, list_cs, list_es), key=lambda pair: abs(
        pair[1][1, 1]) + abs(pair[1][1, 0]))

    # Unzip back into separate sorted lists
    sorted_a, sorted_c, sorted_e = zip(*sorted_pairs)

    # Convert tuples to lists
    sorted_a = list(sorted_a)
    sorted_c = list(sorted_c)
    sorted_e = list(sorted_e)

    for i in range(len(sorted_a)):
        # get dimensions of section
        vecs, x_vals, m0, m1, x0, y0, dy_x, dy_x_frac = setup(
            sorted_a[i], sorted_c[i], sorted_e[i], vecs0, dx, dx_frac)
        print("i = " + str(i), "j = " + str(j))

        if dy == -1:
            dy = dy_x
            dy_frac = dy_x_frac

        # create a dataframe with winning vector at certain points in the section
        df = winners(vecs, x_vals, m0, m1, y0, dx, dy, dx_frac, dy_frac)
        # plot poincare section and save
        try:
            plot(df, vecs, sorted_c[i], j, folder, test=False)
        except Exception as error:
            print(error)
            continue
        df.to_csv(os.path.join(
            "results", folder, f"df_{j}.csv"), index=False)

        # make section object that define winning vector and line equations for boundaries of subsections
        if estimated_output:
            sec_list_estimated = sec_setup_estimated(df, dy)
            secs_estimated = sec_comp_estimated(sec_list_estimated, dx)
            
            with open(os.path.join("results", folder, "secs_estimated_" + str(j) + ".dill"), 'wb') as f:
                dill.dump(secs_estimated, f)
                
            times = time_comp(secs_estimated)
            # plot the pdf for each cusp
            pdf(list(df["time"]), times, dx*2, folder, j)
                
        sec_list_integrals, vec_dict = sec_setup_integrals(df, dy)
        secs_integrals = sec_comp_integrals(df, sec_list_integrals, vec_dict, dx, m1, y0)
        
        with open(os.path.join("results", folder, "secs_integrals_" + str(j) + ".dill"), 'wb') as f:
            dill.dump(secs_integrals, f)

        print(f"section {j} done")
        break


# Code from Sunrose
# Gives non-visibility tori for testing
D = OrigamiDatabase()
q = D.query()
qlist = q.list()


def unit_hor_saddle(O):
    count = 0
    for vert in O.vertices():
        tup = vert.up_right_tuple()
        for i in tup:
            for vert2 in O.vertices():
                tup2 = vert2.up_right_tuple()
                if O.r()(i) in tup2:
                    return True
    return False


def is_unobstructed(O):
    cusp_reps = O.teichmueller_curve().cusp_representatives()
    for item in cusp_reps:
        if not unit_hor_saddle(item[0]):
            return False
    return True


def obstructed(n, **kwargs):
    obstructed = []
    count_obstructed = 0
    p = D.query(nb_squares=n, **kwargs)
    for item in p:
        if not is_unobstructed(item):
            obstructed.append(item)
            count_obstructed += item.teichmueller_curve().orbit_graph().num_verts()
    return (obstructed, count_obstructed)

# list of permutations


def perms_list(n, **kwargs):
    obstructed = []
    p = D.query(nb_squares=n, **kwargs)
    for item in p:
        if not is_unobstructed(item):
            obstructed.append(item)
            for perm in item.teichmueller_curve():
                obstructed.append(perm)
    return obstructed

# get number of new processes to be made when running script_winners


def pool_num(len_alphas):
    num_cores = os.cpu_count()
    num_pools = min(int(num_cores*.5), len_alphas)
    num_loops = len_alphas // num_pools
    print(f"pools: {num_pools}, loops: {num_loops}")
    return num_pools, num_loops

# find the generators of each cusp of the STS


def generators(perm, vecs0):
    generators = []
    a = perm.veech_group().cusps()
    for item in a:
        m = perm.veech_group().cusp_data(item)[0]
        generators.append(m.matrix())
    return generators


class DetailsError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


# For each cusp of the square tiled surface, computes a generating matrix and
# its eigenvectors, then finds the shortest saddle connection in the direction # of the eigenvector by producing a big list of saddle connections and then checking all of them.
# Naming convention: C0 is C^-1, S0 is S^-1
# inputs:
#   perm: a permutation that defines a square tiled surface
#   vecs0: list of saddle connections on this STS to check (this is large)
#
# outputs:
#   alphas: alpha_values
#   Cs: C-matricies
#   C_invs: inverse C-matricies
#   eigs: eigenvalues
#   Ms: C*generator*C^inv
#   generators: generators
#   eigenvecs: eigenvectors

def poincare_setup(perm, vecs0, generators):
    M = mathematica
    eigs = []
    eigenvecs = []
    # get components of each generator
    for matrix in generators:
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]

        # make matrix string for use in mathematica
        matrix_string = f"A = {{{{{a}, {b}}}, {{{c}, {d}}}}}"
        A = M(matrix_string)

        # get eigenvalues using mathematica to get higher precision
        eigenvalues = M('eigenvalues = Eigenvalues[A]')
        eig1 = int(M('eigenvalues[[1]]'))
        eig2 = int(M('eigenvalues[[2]]'))

        # mkae sure they are 1
        if eig1 == eig2:
            if eig1 == 1:
                eigs.append(eig1)
            else:
                raise ValueError("Eigenvalue not equal to 1")
        else:
            raise ValueError("Different eigenvalues")

        # get eigenvectors
        eigenvectors = M('eigenvectors = Eigenvectors[A]')
        # make all the entries integers for easier scaling and grab the first one (matrix isnt diagonlaizable so there is only one unique vector)
        scaledEigenvector = M(
            'scaledEigenvectors = Map[LCM @@ Denominator[#]*# &, eigenvectors][[1]]')
        # get x and y-coords
        scaled_x = int(scaledEigenvector.sage()[0])
        scaled_y = int(scaledEigenvector.sage()[1])
        # create vector
        vec = np.array([[scaled_x], [scaled_y]])
        eigenvecs.append(vec)

    # use numpy to parallelize calculations over vecs0 which is a large list of vectors (around 2.3 million)
    def get_magnitude_slope_sign(vectors):
        # Stack the list of 2D numpy arrays into a 2D array
        vectors = np.hstack(vectors)  # Concatenate into a 2D array

        magnitudes = np.linalg.norm(vectors, axis=0)

        # Make sure the output array for division is explicitly float
        slopes = np.divide(vectors[1], vectors[0], out=np.full_like(
            vectors[1], float('inf'), dtype=np.float64), where=vectors[0] != 0)

        x_signs = np.sign(vectors[0])
        y_signs = np.sign(vectors[1])

        return magnitudes, slopes, x_signs, y_signs

    # get the slope, sign and mag info for the eigenvectors
    eigen_mags, eigen_slopes, eigen_x_signs, eigen_y_signs = get_magnitude_slope_sign(
        eigenvecs)
    # get the slope, sign and mag info for all saddle connections
    saddle_mags, saddle_slopes, saddle_x_signs, saddle_y_signs = get_magnitude_slope_sign(
        vecs0)

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

        valid_saddles = np.where(
            slope_matches & x_sign_matches & y_sign_matches)[0]

        if len(valid_saddles) == 0:
            raise DetailsError(
                f"No saddle connection found for eigenvector {eigenvecs[i]} in provided vectors data")

        # Select saddle with smallest magnitude
        smallest_idx = valid_saddles[np.argmin(saddle_mags[valid_saddles])]
        saddle_vecs.append(vecs0[smallest_idx])

    saddle_vecs = np.array(saddle_vecs)

    # find c_inv and c
    # C-matricies
    Cs = []
    # C-inverses
    C0s = []
    # alpha values for sections
    alphas = []
    # jordan decomposition with pre-scaled matricies
    Js = []
    # jordan decomposition with scaled matricies, used for alphas
    Ms = []
    # pre-scaled C-matrix
    Ss = []
    # pre-scaled C^-1 matrix
    S0s = []
    # factor needed to take eigenvec -> saddle connection
    mults = []
    # 2x2 matrix using mults to take S0 and S to C0 and C
    scales = []
    for matrix, saddle_vec, eigen_vec in zip(generators, saddle_vecs, eigenvecs):
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]

        matrix_string = f"A = {{{{{a}, {b}}}, {{{c}, {d}}}}}"
        A = M(matrix_string)

        # S0 is pre-scaled C^-1 matrix, S*G*S0 = J
        S0, J = M('{S0, J} = JordanDecomposition[A]')
        S = M('S = Inverse[S0]')

        # print(S0)
        # print(S)
        # print()
        # print("-----------------------")

        # convert S to a usable version in python
        S_np = np.array(S.sage())
        a_ = frac(str(S_np[0][0]))
        b_ = frac(str(S_np[0][1]))
        c_ = frac(str(S_np[1][0]))
        d_ = frac(str(S_np[1][1]))
        S = np.array([[a_, b_], [c_, d_]], dtype=object)

        # convert S0 to a usable version in python
        S0_np = np.array(S0.sage())
        a_2 = frac(str(S0_np[0][0]))
        b_2 = frac(str(S0_np[0][1]))
        c_2 = frac(str(S0_np[1][0]))
        d_2 = frac(str(S0_np[1][1]))
        S0 = np.array([[a_2, b_2], [c_2, d_2]], dtype=object)

        # convert J to a usable version in python
        M_np = np.array(J.sage())
        a_1 = int(M_np[0][0])
        b_1 = int(M_np[0][1])
        c_1 = int(M_np[1][0])
        d_1 = int(M_np[1][1])
        J = np.array([[a_1, b_1], [c_1, d_1]])

        # make sure J is the correct Jordan decomposition
        if a_1 != 1 or c_1 != 0 or d_1 != 1:
            raise ValueError(f"wrong J: {J}")

        # get determinant of S0
        detC = (a_2*d_2) - (b_2*c_2)
        factor = int(1)/sqrt(detC)
        # normalize S0, call is c0
        c0 = (factor * S0)

        # find factor needed to send (1,0) to the saddle connection
        x_ = (c0@np.array([[1], [0]]))[0][0]
        y_ = (c0@np.array([[1], [0]]))[1][0]
        if x_ != 0:
            mult = saddle_vec[0][0]/x_
        elif y_ != 0:
            mult = saddle_vec[1][0]/y_
        else:
            raise ValueError(
                "Both coordinates in saddle_vec are zero. Division by zero is not defined.")

        # create a determinant 1 matrix that will scale c0 to send (1,0) to the saddle connection with correct length
        scale = np.array([[mult, 0], [0, int(1)/mult]], dtype=object)
        c0 = c0 @ scale
        # find the inverse
        c = np.array(
            [[c0[1][1], -c0[0][1]], [-c0[1][0], c0[0][0]]], dtype=object)

        # compute new Jordan decomposition
        J_new = c@matrix@c0
        a_1 = int(J_new[0][0])
        b_1 = int(J_new[0][1])
        c_1 = int(J_new[1][0])
        d_1 = int(J_new[1][1])

        # ensure its correctly formatted
        if a_1 != 1 or c_1 != 0 or d_1 != 1:
            raise ValueError(f"wrong J_new: {J_new}")

        Cs.append(c)
        C0s.append(c0)
        Ms.append(J_new)
        Js.append(J)
        S0s.append(S0)
        Ss.append(S)
        scales.append(scale)
        mults.append(mult)
        alphas.append(J_new[0][1])

    return alphas, Cs, C0s, eigs, Ms, generators, eigenvecs

# A try-catch wrapper on poincare_details, for testing purposes.


def try_poincare_setup(sts_data, trys):
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

# compute the dimensions of the Poincare section
# using the section vector guarantees a finite number of winners: https://arxiv.org/abs/2102.10069v2

# input:
    # alpha: alpha value for the section
    # c: c-matrix for the section
    # eig: eigenvalue for the section
    # vecs0: original set up saddle connections for the STS
    # dx: x-spacing used for point sampling in the function "winners"

# output:
    # vecs: new vectors to be used for computations, acted on by the c-matrix for the given cusp
    # x_vals: x-values used for computation
    # m0: slope for the top portion of the section
    # m1: slope for the bottom portion of the section
    # x0: x component of the "section" vector
    # y0: y component of the "section" vector
    # dx_y: y-spacing used for point sampling in the function "winners"


def setup(alpha, c, eig, vecs0, dx, dx_frac):
    x_vals = np.arange(dx, 1, dx)
    # create list of new vectors
    vecs1 = c@vecs0
    vecs = []
    for item in vecs1:
        vecs.append(item)

    # section vector is the vector with smallest y-component with y > 0 and corresponding shortest x-component where x > 0. It defines the y-intercept of the section and partially defines slopes of the lines of the section
    sec_vec = np.array([[None], [None]])
    for i in range(len(vecs)):
        # check if saddle connection has negative components
        if (vecs[i][0][0] <= 0 or vecs[i][1][0] <= 0):
            continue
        else:
            # establish a section vector if there is none
            if sec_vec[0][0] == None:
                sec_vec = vecs[i]
            # compare section vector to saddle connection
            else:
                # if sec_vec y-component is larger, define section vector as this saddle connection
                if sec_vec[1][0] > vecs[i][1][0]:
                    sec_vec = vecs[i]
                # if y-components are the same, compare x-components
                elif sec_vec[1][0] == vecs[i][1][0]:
                    if sec_vec[0][0] > vecs[i][0][0]:
                        sec_vec = vecs[i]
    # if no such section vector exists, assign the vector (0,1) which is the same as not having a section vector in calculation
    if sec_vec[0][0] == None:
        sec_vec = np.array([[0], [1]])
    # assign components to variables
    x0 = sec_vec[0][0]
    y0 = sec_vec[1][0]

    # slopes of top and bottom lines of section
    m0 = -x0/y0
    m1 = -(x0/y0 + alpha)

    # defines vertical step when finding winners
    dy_x = (m0 - m1) * dx
    dy_x_frac = frac((m0 - m1) * dx_frac)
    
    return vecs, x_vals, m0, m1, x0, y0, dy_x, dy_x_frac


# It first computes winning vectors from all possible vectors on points on the
# verticals of the poincare section. Then it computes winners
# for all points in the poincare section from the subsetset of winners calculated in the first step.
#
# inputs:
#   vecs: array of vectors, after being multiplied by C matrix
#   x_vals: array of x-values in the poincare section, all between [0,1]
#   m0: slope of the top line of the poincare section
#   m1: slope of the bottom line of the poincare section
#   y0: 1/y0 is the y-coordinate of the top left corner point of the section
#   dx: step size for separation of points in the x-direction to sample in the section
#   dy: step size for separation of points in the y-direction to sample in the section
#
# output:
    # df: dateframe that has the x and y coordinates that were sampled, the winning saddle connection at those coordinates, the label given to the vector for plotting purposes, and the rounded return time associated with that point

def winners(vecs0, x_vals, m0, m1, y0, dx, dy, dx_frac, dy_frac):
    # dictionary for plotting
    saddle_dict = {}
    saddle_dict["x"] = []
    saddle_dict["y"] = []
    saddle_dict["lab"] = []
    saddle_dict["vec"] = []
    saddle_dict["time"] = []
    possible_vecs = []
    winners = []
    vecs_float = np.hstack([arr.astype(float) for arr in vecs0])  # float operations
    vecs_frac = np.hstack(vecs0)  # fraction operations
    
    # verticals
    t0 = time()
    x_vals_float = [0.9, 0.95, 0.9995]
    x_vals_frac = [frac('90/100'), frac('95/100'), frac('9995/10000')]
    
    y_vals_float_total = []
    y_vals_frac_total = []
    for x_float, x_frac in zip(x_vals_float, x_vals_frac):
        y_vals_float = []
        y_vals_frac = []
        current_y = m1*x_frac + frac(int(1), y0) + dy_frac
        while current_y <= m0*x_frac + frac(int(1), y0) - dy_frac:
            y_vals_frac.append(current_y)
            y_vals_float.append(float(current_y))
            current_y = current_y + dy_frac
        y_vals_frac_total.append(y_vals_frac)
        y_vals_float_total.append(y_vals_float)
    
    print(f"number of verticals: {len(x_vals_float)}")
    for i, a in enumerate(x_vals_float):
        for j, b in enumerate(y_vals_float_total[i]):
            Mab_float = np.array([[a, b], [0, 1/a]], dtype='float')
            # Apply the transformation to all vectors at once
            new_vecs = Mab_float @ vecs_float
    
            new_vecs_x = new_vecs[0, :]
            new_vecs_y = new_vecs[1, :]
    
            # Filter based on conditions (x > 0, y/x > 0, and x <= 1)
            filtered_indices = (new_vecs_x > 0) & (new_vecs_x <= 1) & (
                new_vecs_y / np.where(new_vecs_x == 0, np.inf, new_vecs_x) > 0)
    
            # get x and y components of filtered vectors
            filtered_x = new_vecs_x[filtered_indices]
            filtered_y = new_vecs_y[filtered_indices]
            # Filter the fractional vecs also
            filtered_vecs_frac = vecs_frac[:, filtered_indices]
    
            if filtered_x.size == 0:
                winners.append(None)
                continue
    
            # Calculate slopes
            slopes = filtered_y / filtered_x
    
            # numerical smallest slope, could have error
            min_slope = np.min(slopes)
            # find (the indices of vectors with) similar numeric slopes in casae something got messed up
            nearby_slopes = np.isclose(slopes, min_slope)
            nearby_slope_indices = np.where(nearby_slopes)[0]  # Get the integer indices
            # go get those indices from the filtered vectors with frac representation
            vecs_frac_close = filtered_vecs_frac[:, nearby_slope_indices]       
    
            Mab_frac = np.array([[x_vals_frac[i], y_vals_frac_total[i][j]], [frac(int(0)), frac(int(1), x_vals_frac[i])]], dtype='object')
            candidates_frac = Mab_frac @ vecs_frac_close
            slopes_frac = candidates_frac[1, :] / candidates_frac[0, :]
    
            winner_idx = None
            for idx, slope, x_comp in zip(nearby_slope_indices, slopes_frac, candidates_frac[0, :]):
                if winner_idx == None:
                    winner_idx = idx
                    current_slope = slope
                    current_x_comp = x_comp
                else:
                    if slope > current_slope:
                        continue
                    elif slope < current_slope:
                        winner_idx = idx
                        current_slope = slope
                        current_x_comp = x_comp
                    else:
                        if x_comp < current_x_comp:
                            winner_idx = idx
                            current_slope = slope
                            current_x_comp = x_comp
                        else:
                            continue
            # Among candidates with minimal slope, pick the one with smallest x
            winner = filtered_vecs_frac[:, winner_idx]
            winners.append(winner.reshape(2, 1))
    t1 = time()
    print("verticals done: " + str(t1 - t0))
    # Step 1: Filter out None values efficiently
    winners_filtered = [w for w in winners if w is not None]
    
    # Step 2: Find unique vectors using a dictionary (faster than list searching)
    unique_dict = defaultdict(list)
    for winner in winners_filtered:
        # Convert to tuple for hashability
        key = tuple(winner.flatten())
        unique_dict[key].append(winner)
    
    # Step 3: Get first occurrence of each unique vector
    possible_vecs = [vectors[0] for vectors in unique_dict.values()]

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
    possible_vecs_float = [vec.astype(float) for vec in possible_vecs]

    # print(possible_vecs)
    # print(possible_vecs_float)

    for i in range(len(possible_vecs_float)):
        # apply Mab matrix, perform horocycle flow and find time t to horizontal
        a = horo@(Mab@possible_vecs_float[i])
        # print(a)
        t_dict[i] = lambdify([x, y], solve(a[1][0], t)[0])

    # for each point (a,b) in the poincare section, apply the Mab matrix to each vector and look for "winners". Winners have smallest possible slope that is greater than zero and 0 < x-component <= 1
    for a in x_vals:
        y_vals = np.arange(m1*a + 1/y0 + dy, m0*a + 1/y0 - dy, dy)
        for b in y_vals:
            check = 0
            winner_slope = None
            winner = None
            Mab = np.array([[a, b], [0, 1/a]], dtype='float')
            for vec, vec1 in zip(possible_vecs, possible_vecs_float):
                new = Mab@vec1
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
                    # this is an error
                    if saddle_dict["time"][-1] < 0:
                        saddle_dict["time"][-1] = 1000
            # if there is no winner, at (a,b), add a label so the df and plot can still be made. These section will later be made blank for trouble-shoooting
            if check == 0:
                saddle_dict["lab"].append(len(vecs0))
                saddle_dict["time"].append(1000)

    df = pd.DataFrame.from_dict(saddle_dict)
    return df

# plot the Poincare sections and save them
def plot(df, vecs, c, j, folder, test=False):
    fig, ax = plt.subplots(figsize=(10, 10))
    # plot winners
    ax.scatter(df[df["lab"] != len(vecs)]["x"], df[df["lab"] != len(
        vecs)]["y"], c=df[df["lab"] != len(vecs)]["lab"], s=0.1)
    # for points with no winner, make white
    ax.scatter(df[df["lab"] == len(vecs)]["x"], df[df["lab"]
               == len(vecs)]["y"], c="white", alpha=0.5, s=0.1)
    ax.set_aspect('auto', adjustable='box')

    # plot outline of poincare section
    min_x = min(df["x"])
    y1 = max(df[df["x"] == min_x]["y"])
    y2 = min(df[df["x"] == min_x]["y"])
    max_x = max(df["x"])
    y3 = max(df[df["x"] == max_x]["y"])
    y4 = min(df[df["x"] == max_x]["y"])
    ax.plot([min_x, max_x, max_x, min_x, min_x],
            [y1, y3, y4, y2, y1], c="black")

    plt.title(str(c) + "\n", loc="left", fontsize=25)
    if test == False:
        # take out
        plt.show()
        plt.savefig(os.path.join(
            "results", folder, "section_" + str(j)))
    # for troubleshooting
    if test == True:
        # display vectors on the right edge of the section from top to bottom
        vecs_list = list(df[df["x"] == max(df["x"])]["vec"])
        unique_arrays = set(tuple(map(tuple, arr)) for arr in vecs_list)
        unique_arrays = [np.array(arr) for arr in unique_arrays]
        unique_arrays.reverse()
        for arr in unique_arrays:
            print(tuple(item[0] for item in arr))
    # if there are points with no winning vec, raise an error
    if len(df[df["lab"] == len(vecs)]) != 0:
        raise ValueError("Poincare section has empty portion")
    plt.show()
    plt.close(fig)

class Section:
    def __init__(self):
        self.vec = None
        # equations
        self.top = []
        self.bottom = []
        # points
        self.points_bottom = None

# this is used for computing the integral-based pdfs for a given poincare section
def sec_setup_integrals(df, dy):
    sec_list = []
    labs = df["lab"].unique()
    # dictionary for label of vector to the vector itself
    vec_dict = {}
    # indexing through the winning saddle connections that appear in the section
    for lab in labs:
        sec_dict = {}
        # get the subsection defined by one winner
        df1 = df[df["lab"] == lab]
        # add vector to dictionary
        vec_dict[lab] = df1['vec'].iloc[int(0)]
        # get all x values for the subsection
        xs = df1["x"]
        # sort them into a unique list
        xs = sorted(list(set(xs.tolist())))
        # list of the y-values that make up the bottom portion of the subsection
        y_bottoms = []
        for x in xs:
            # for a given "x" find the top and bottom y-values
            y_top = max(df1[df1["x"] == x]["y"])
            y_bottom = min(df1[df1["x"] == x]["y"])
            # ensures the section is convex, not concave
            if len(df1[df1["x"] == x]["y"]) < (y_top - y_bottom)/dy:
                print("len: " + str(len(df1[df1["x"] == x]["y"])))
                print("ytop: " + str(y_top))
                print("ybottom: " + str(y_bottom))
                print("dy: " + str(dy))
                print(x)
                print(df1[df1["x"] == x]["y"])
                raise ValueError(
                    "Section has more than 2 points for a given 'x'")
            y_bottoms.append(y_bottom)
        sec_dict['x'] = xs
        sec_dict['bottom'] = y_bottoms
        sec_dict['vec'] = df1['vec'].iloc[int(0)]
        sec_list.append(sec_dict)
    return sec_list, vec_dict

# this is used for computing the integral-based pdfs for a given poincare section
def sec_comp_integrals(df, sec_list, vec_dict, dx, m1, y0):
    secs = []
    problem_xs = []
    # go through each section and find its left-most point and append it to a list
    for i in range(len(sec_list)):
        xs = sec_list[i]['x']
        problem_xs.append(xs[0])
    problem_xs = sorted(list(set(problem_xs)))
    # go through all these problem_xs and get the unique values. Due to sampling, there mayb be two subsections that should begin at x = 0 but begin at x = 0.01 and x = 0.0105. Count these points as the same and only add one of them
    problems = [problem_xs[0]]
    for i in range(len(problem_xs) - 1):
        # these points are the same but may not be equal due to sampling methods and python rounding
        if abs(problem_xs[i] - problem_xs[i+1]) <= 0.01:
            continue
        problems.append(problem_xs[i+1])
    problems.append(1)
    
    # find all the midpoints between these problem points
    prob_mids = []
    for i in range(len(problems) - 1):
        prob_mids.append((problems[i] + problems[i+1])/2)

    # go through each midpoint and find the closest x-value that was sampled in our winners df. 
    i = 0
    use_points = []
    for item in np.arange(dx, 1, dx):
        if item < prob_mids[i]:
            continue
        i += 1
        use_points.append(item)
        if i == len(prob_mids):
            break

    for i in range(len(sec_list)):
        xs = sec_list[i]['x']
        bottom = sec_list[i]['bottom']
        sec = Section()
        sec.vec = sec_list[i]['vec']

        # list to accumulate what sections are directly about the given subsection. These will be used to determine the equations of the lines that make up the bottom part of the subsection
        bottom_lab_list = []
        for x_, y_ in zip(xs, bottom):
            if x_ not in use_points:
                continue
            # list of y-values at the given x-value
            possible_ys = sorted(list(df[(df["x"] == x_)]["y"]), reverse=True)
            # find the index of the y-value that is on the bottom edge of the current subsection at the current x-value
            y_index = possible_ys.index(y_)

            # Check if there is a y-value in our winners df that is below the y_index value. If there is, find out what the label/vector is in that different subsection
            if y_index + 1 < len(possible_ys):
                # Get the value at index + 1
                next_y = possible_ys[y_index + 1]
                filtered_df = df[(df["x"] == x_) & (df["y"] == next_y)]
                # Get the value in the "lab" column of the first row
                lab_value = filtered_df["lab"].iloc[0]
                if lab_value not in bottom_lab_list:
                    bottom_lab_list.append(lab_value)
            # if there are no y-values, then this subsection has a bottom edge on the boundary of the entire poincare section. Record this label as -1
            else:
                if -1 not in bottom_lab_list:
                    bottom_lab_list.append(int(-1))

        # list of equations for the bottom segments of our subsection. The equations of the top segemnts of a given subsection is determined by the winning saddle connection. Use this fact to get the equations of the bottom segments that share an edge with these top segments
        eqs = []
        for lab in bottom_lab_list:
            if lab != -1:
                vec = vec_dict[lab]
                # for a winning saddle connection (x0,y0), the equation of its top line is (-x0/y0)*x + 1/y0
                eqs.append(-frac(vec[0][0], vec[1][0])
                           * x + frac(int(1)/vec[1][0]))
            # if there is no other subsection below the current subsection, the equation for the bottom segment is the boundary of the poincare section
            else:
                eqs.append(frac(m1) * x + frac(int(1)/y0))

        # get the exact x-values for where the equation of the bottom segments change
        # points where the bottom change
        points = [xs[0]]
        for i in range(len(eqs) - 1):
            eq1 = eqs[i]
            eq2 = eqs[i + 1]
            # Set eq1 = eq2 and solve for x
            intersection_x = solve(Eq(eq1, eq2), x)
            if intersection_x:  # If there is a solution
                points.append(intersection_x[0])
            else:
                raise ValueError(f"intersection between bottom equation {i} and bottom equation {i+1} point failed")

        # define a top function for the top of the subsection
        sec.top.append(-frac(sec.vec[0][0], sec.vec[1][0])
                       * x + frac(int(1)/sec.vec[1][0]))

        # modifying the first entry in points so that it is not the rounded version. Finding the x-component of the intersection between the "top" equation and the first "bottom" equation
        eq1 = sec.top[0]
        eq2 = eqs[0]

        intersection_x = solve(Eq(eq1, eq2), x)
        if intersection_x:  # If there is a solution
            points[0] = intersection_x[0]
        else:
            raise ValueError("intersection between top and first bottom point failed")

        for eq in eqs:
            sec.bottom.append(eq)
            sec.points_bottom = points
        secs.append(sec)
    return secs

def save_arrays_to_file(file_path, arrays_list):
    # Save arrays to a single NumPy file
    np.save(file_path, arrays_list)


def load_arrays_from_file(file_path):
    # Load arrays from the NumPy file
    arrays_list = np.load(file_path, allow_pickle=True)
    # Ensure each element in the list is a NumPy array
    arrays_list = [np.array(array) for array in arrays_list]
    return arrays_list

def read_df(folder, cusp):
    df = pd.read_csv(os.path.join(
        "results", folder, "df_" + str(cusp)))

    def read_vec(s):
        s = s.replace("[", "")
        s = s.replace("]", "")
        nums = s.split("\n")
        nums[0] = frac(nums[0])
        nums[1] = frac(nums[1])
        return np.array([[nums[0]], [nums[1]]])

    df["vec"] = df["vec"].apply(read_vec)
    return df

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class ComputationsTestSuite(unittest.TestCase):
    """Basic test cases."""

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
        t0 = time()
        details_1 = try_poincare_details((perm, vecs), 3)
        t1 = time()
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
        t0 = time()
        details_1 = try_poincare_details((perm, vecs), 1)
        t1 = time()
        print(f'  runtime: {t1-t0}')
        # print(details_1)
        # output = compute_on_cusp(details, vecs)

    def test_compute_winning_vecs_on_edges(self):
        # this is permutations[3] from generate_permutations(7)
        perm = Origami(
            '(1)(2)(3)(4,5)(6,7)', '(1,2,3,4,6)(5,7)')
        # this is corresponding vector data computed with full library
        vectors = load_arrays_from_file(os.path.join(
            "Poincare_Sections", "vecs", "vecs7-3.npy"))[0:1000]

        generators = []
        a = perm.veech_group().cusps()
        # TODO: fix a  to actually have validatable output data = [Infinity, -4/7, -1/5, -5/29, -3/5, -15/22, 0, -2/3, -1/2, -1/6]
        print(f'cusps:{a}')
        # print(type(a[0]))
        for item in a:
            m = perm.veech_group().cusp_data(item)[0]
            generators.append(m.matrix())

        # mock data for setup
        alphas = [1.0, 1076.0, 104.0, 106.0, 5862.0,
                  663.0, 10.0, 120.0, 455.0, 119.0]
        c_matrices = [np.array([[1., 0.], [0., 1.]]), np.array([[0.00557628, -0.04275121], [0.04275121,  0.00557628]]), np.array([[0.03846305, -0.19231131], [0.19231131,  0.03846305]]), np.array([[0.03773585, -0.13207547], [0.13207547,  0.03773585]]), np.array([[0.00409406, -0.0317294], [0.0317294,  0.00409406]]),
                      np.array([[0.02262495, -0.06334904], [0.06334904,  0.02262495]]), np.array([[0.,  0.5], [-0.5, -0.]]), np.array([[0.1000019, -0.30000253], [0.30000253,  0.1000019]]), np.array([[0.01538339, -0.12307211], [0.12307211,  0.01538339]]), np.array([[0.05882581, -0.2352984], [0.2352984,  0.05882581]])]
        eigenvectors = [matrix([[1, 1], [0, 1]]), matrix([[139, 18], [-1058, -137]]), matrix([[21, 4], [-100, -19]]), matrix([[29, 8], [-98, -27]]), matrix([[745,    96], [-5766,  -743]]), matrix(
            [[211, 75], [-588, -209]]), matrix([[1,   0], [-10,   1]]), matrix([[37,   12], [-108,  -35]]), matrix([[57,    7], [-448,  -55]]), matrix([[29,    7], [-112,  -27]])]
        # last_piece_of_details = [np.array([[1], [0]]), np.array([[1.], [-7.66666667]]), np.array([[1], [-5]]), np.array([[1.], [-3.5]]), np.array([[1.], [-7.75]]),
        #      np.array([[1.], [-2.8]]), np.array([[0], [1]]), np.array([[1], [-3]]), np.array([[1], [-8]]), np.array([[1], [-4]])]

        n_squares = 7
        dx = 0.0005
        index = 4

        expected_len = [1, 3, 1, 1, 3, 1, 44, 1, 1, 1]
        # TODO: test against expected values,
        # expected_winners = [[np.array([[0], [2]])],
        #                     [np.array([[-6.76584374], [0.79554982]]), np.array([[-5.56323358], [0.68216519]]), np.array([[-0.08550242], [0.01115256]])]]
        for j in range(0, 10):
            vecs, x_vals, m0, m1, x0, y0, dx_y = setup(
                alphas[j], c_matrices[j], eigenvectors[j], vectors, dx, False)

            t0 = time()
            result = winners(
                vecs, x_vals, m0, m1, y0, dx, dx_y)
            t1 = time()
            print(f'time: {t1-t0}')
            print(f'edge winning vecs: {result}')
            print(f'len: {len(result)}')
            self.assertEqual(len(result), expected_len[j])

            # numpy doesn't like mock data :(
            # for i in range(0, len(result)):
            #     self.assertTrue(np.array_equal(
            #         result[i], expected_winners[j][i]))
