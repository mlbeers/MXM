from sage.all import *
from sage.all import matrix  # testing
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
#from flatsurf import *
import pwlf
import os
#import math
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

# run_script - brief description
# inputs:
# vecs0 - array of ... numpy
# a - the alpha value... 
def compute_poincare_sections(vecs0, a, c, e, dx, j, n_squares, index):
    print(f"id: {os.getpid()}")
    list_as = []
    list_cs = []
    list_es = []
    for a_, c_, e_ in zip(a, c, e):
        list_as.append(a_[j])
        list_cs.append(c_[j])
        list_es.append(e_[j])
    
    # Sort a, c, and e together based on the absolute value of the lower-right entry of c
    sorted_pairs = sorted(zip(list_as, list_cs, list_es), key=lambda pair: abs(pair[1][1, 1]) + abs(pair[1][1,0]))
    
    # Unzip back into separate sorted lists
    sorted_a, sorted_c, sorted_e = zip(*sorted_pairs)
    
    # Convert tuples to lists
    sorted_a = list(sorted_a)
    sorted_c = list(sorted_c)
    sorted_e = list(sorted_e)

    for i in range(len(sorted_a)):
        # get dimensions of section
        vecs, x_vals, m0, m1, x0, y0, dx_y, z = setup(
            sorted_a[i], sorted_c[i], sorted_e[i], vecs0, dx, True)
        print("i = " + str(i), "j = " + str(j))

        if float(z) <= float(1/50000):
            print("too small")
            continue

        # create a dataframe with winning vector at certain points in the section
        df = winners1(vecs, x_vals, m0, m1, y0, dx, dx_y)
        # plot poincare section and save
        try:
            plot(df, vecs, sorted_c[i], j, n_squares, index, test=False)
        except Exception as error:
            print(error)
            continue
        df.to_csv(os.path.join(
            "results", f"{n_squares} - {index}", f"df - {j}.csv"), index=False)

        # make section object that define winning vector and line equations for boundaries of subsections
        sec_list = sec_setup(df, dx_y)
        secs = sec_comp(sec_list, dx)
        sec_list2, vec_order, vec_dict = sec_setup2(df, dx_y)
        secs2 = sec_comp2(df, sec_list2, vec_order, vec_dict, dx, dx_y, m1, y0)

        
        with open(os.path.join("results", f"{n_squares} - {index}", "secs - " + str(j) + ".dill"), 'wb') as f:
            dill.dump(secs, f)

        with open(os.path.join("results", f"{n_squares} - {index}", "secs_integrals - " + str(j) + ".dill"), 'wb') as f:
            dill.dump(secs2, f)

        times = time_comp(secs)

        # plot the pdf for each cusp
        pdf(list(df["time"]), times, dx*2, n_squares, index, j)
        
        print(f"section {j} done")

        break

def pool_num(len_alphas):
    num_cores = os.cpu_count()
    num_pools = min(int(num_cores*.5), len_alphas)
    num_loops = len_alphas // num_pools
    print(f"pools: {num_pools}, loops: {num_loops}")
    return num_pools, num_loops

def generators(perm, vecs0):
    # find the generators of each cusp of the STS
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

def poincare_setup(perm, vecs0, generators):
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
    
    # find c_inv and c
    Cs = []
    C0s = []
    alphas = []
    Ms = []
    Ss = []
    S0s = []
    Ss = []
    scales = []
    mults = []
    Js = []
    for matrix, saddle_vec, eigen_vec in zip(generators, saddle_vecs, eigenvecs):
        a = matrix[0][0]
        b = matrix[0][1]
        c = matrix[1][0]
        d = matrix[1][1]
    
        matrix_string = f"A = {{{{{a}, {b}}}, {{{c}, {d}}}}}"
        A = M(matrix_string)
    
        S0, J = M('{S0, J} = JordanDecomposition[A]')
        S = M('S = Inverse[S0]')

        # print(S0)
        # print(S)
        # print()
        # print("-----------------------")
    
        S_np = np.array(S.sage())
        a_ = frac(str(S_np[0][0]))
        b_ = frac(str(S_np[0][1]))
        c_ = frac(str(S_np[1][0]))
        d_ = frac(str(S_np[1][1]))
        S = np.array([[a_, b_], [c_, d_]], dtype=object)
    
        M_np = np.array(J.sage())
        a_1 = int(M_np[0][0])
        b_1 = int(M_np[0][1])
        c_1 = int(M_np[1][0])
        d_1 = int(M_np[1][1])
        J = np.array([[a_1, b_1], [c_1, d_1]])

        S0_np = np.array(S0.sage())
        a_2 = frac(str(S0_np[0][0]))
        b_2 = frac(str(S0_np[0][1]))
        c_2 = frac(str(S0_np[1][0]))
        d_2 = frac(str(S0_np[1][1]))
        S0 = np.array([[a_2, b_2], [c_2, d_2]], dtype=object)
    
        if a_1 != 1 or c_1 != 0 or d_1 != 1:
            raise ValueError(f"wrong J: {J}")

        detC = (a_2*d_2) - (b_2*c_2)
        factor = int(1)/sqrt(detC)
        c0 = (factor * S0)
        x_ = (c0@np.array([[1],[0]]))[0][0]
        y_ = (c0@np.array([[1],[0]]))[1][0]
        mult0 = saddle_vec[0][0]/x_
        mult1 = saddle_vec[1][0]/y_
        if x_ != 0:
            mult = mult0
        elif y_ != 0:
            mult = mult1
        else:
            raise ValueError("Both coordinates in saddle_vec are zero. Division by zero is not defined.")
        scale = np.array([[mult, 0], [0, int(1)/mult]], dtype=object)
        c0 = c0 @ scale
        c = np.array([[c0[1][1], -c0[0][1]], [-c0[1][0], c0[0][0]]], dtype=object)

        J_new = c@matrix@c0
        a_1 = int(J_new[0][0])
        b_1 = int(J_new[0][1])
        c_1 = int(J_new[1][0])
        d_1 = int(J_new[1][1])

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
def try_poincare_details(sts_data, trys):
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


def setup(alpha, c, eig, vecs0, dx, improved=True):
    x_vals = np.arange(dx, 1, dx)
    # for the poincare section with matrix c, the original saddle connections are acted on by this matrix
    vecs1 = c@vecs0
    # find vectors such that -10 <= x <= 10 and 0 < y <= 10
    vecs = []
    for item in vecs1:
        vecs.append(item)

    # section vector is vector with smallest y-component with y > 0 and corresponding x-component with 0 <= x <= 1 with smallest slope. It defines the y-intercept of the section and partially defines slopes of the lines of the section
    sec_vec = np.array([[None], [None]])
    for i in range(len(vecs)):
        if (vecs[i][0][0] <= 0 or vecs[i][1][0] <= 0):
            continue
        else:
            if sec_vec[0][0] == None:
                sec_vec = vecs[i]
            else:
                if sec_vec[1][0] > vecs[i][1][0]:
                    sec_vec = vecs[i]
                elif sec_vec[1][0] == vecs[i][1][0]:
                    if sec_vec[0][0] > vecs[i][0][0]:
                        sec_vec = vecs[i]
    if sec_vec[0][0] == None:
        sec_vec = np.array([[0], [1]])
    x0 = sec_vec[0][0]
    y0 = sec_vec[1][0]
    z = y0

    # slopes of top and bottom lines of section
    if improved == True:
        m0 = -x0/y0
        m1 = -(x0/y0 + alpha)
    else:
        m0 = 0
        m1 = -alpha
        x0 = 0
        y0 = 1

    # defines vertical step when finding winners
    global dx_y
    dx_y = (m0 - m1) * dx

    return vecs, x_vals, m0, m1, x0, y0, dx_y, z


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
def winners(vecs, x_vals, m0, m1, y0, dx, dx_y):
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


def plot(df, vecs, c, j, n_squares, index, test=False):
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
        #take out
        plt.show()
        plt.savefig(os.path.join(
            "results", f"{n_squares} - {index}", "section - " + str(j)))
    # for troubleshooting
    if test == True:
        # display vectors on the right edge of the section from top to bottom
        vecs_list = list(df[df["x"] == max(df["x"])]["vec"])
        unique_arrays = set(tuple(map(tuple, arr)) for arr in vecs_list)
        unique_arrays = [np.array(arr) for arr in unique_arrays]
        unique_arrays.reverse()
        for arr in unique_arrays:
            print(tuple(item[0] for item in arr))
    if len(df[df["lab"] == len(vecs)]) != 0:
        raise ValueError("Poincare section has empty portion")
    plt.show()
    plt.close(fig)


class Section:
    def __init__(self, x, top, bottom):
        self.vec = None
        self.pwlf_top = pwlf.PiecewiseLinFit(x, top)
        self.pwlf_bottom = pwlf.PiecewiseLinFit(x, bottom)
        # equations
        self.top = []
        self.bottom = []
        # lambdified equations
        self.f_top = []
        self.f_bottom = []
        # points
        self.points_top = None
        self.points_bottom = None

    def t(self):
        x, y, t = sym.symbols('x y t')
        Mab = np.array([[x, y], [0, 1/x]])
        horo = np.array([[1, 0], [-t, 1]])
        a = horo@(Mab@self.vec)
        return solve(a[1][0], t)[0]

    # find the time equation in terms of y
    def y(self):
        x, y, t = sym.symbols('x y t')
        Mab = np.array([[x, y], [0, 1/x]])
        horo = np.array([[1, 0], [-t, 1]])
        a = horo@(Mab@self.vec)
        return solve(a[1][0], y)[0]


def sec_setup(df, dx_y):
    sec_list = []
    global labs
    labs = df["lab"].unique()
    for lab in labs:
        sec_dict = {}
        df1 = df[df["lab"] == lab]
        xs = df1["x"]
        xs = sorted(list(set(xs.tolist())))
        y_tops = []
        y_bottoms = []
        for x in xs:
            # for a given "x" find the max and minimum y-values
            y_top = max(df1[df1["x"] == x]["y"])
            y_bottom = min(df1[df1["x"] == x]["y"])
            # ensures the section is convex, not concave
            if len(df1[df1["x"] == x]["y"]) < (y_top - y_bottom)/dx_y:
                print("len: " + str(len(df1[df1["x"] == x]["y"])))
                print("ytop: " + str(y_top))
                print("ybottom: " + str(y_bottom))
                print("dx_y: " + str(dx_y))
                print(x)
                print(df1[df1["x"] == x]["y"])
                raise ValueError(
                    "Section has more than 2 points for a given 'x'")
            y_tops.append(y_top)
            y_bottoms.append(y_bottom)
        y_tops = np.array(y_tops, dtype='float32')
        y_bottoms = np.array(y_bottoms, dtype='float32')
        xs = np.array(xs, dtype='float32')
        sec_dict['x'] = xs
        sec_dict['top'] = y_tops
        sec_dict['bottom'] = y_bottoms
        sec_list.append(sec_dict)
    return sec_list


def sec_comp(sec_list, dx):
    secs = []
    for i in range(len(sec_list)):
        x = sec_list[i]['x']
        top = sec_list[i]['top']
        bottom = sec_list[i]['bottom']
        sec = Section(x, top, bottom)
        sec.vec = label_dict[labs[i]]

        # use piece-wise linear regression to find the equations of the lines for subsection
        # top
        num = 1
        check = True
        while check:
            breaks1 = sec.pwlf_top.fit(num)
            score = sec.pwlf_top.r_squared()
            if score == float("-inf") or 1 - score < dx*2:
                check = False
            if num > 3:
                breaks1 = sec.pwlf_top.fit(1)
                score = sec.pwlf_top.r_squared()
                check = False
            num += 1

        # bottom
        num = 1
        check = True
        while check:
            breaks2 = sec.pwlf_bottom.fit(num)
            score = sec.pwlf_bottom.r_squared()
            if score == float("-inf") or 1 - score < dx*2:
                check = False
            if num > 3:
                breaks2 = sec.pwlf_bottom.fit(num)
                score = sec.pwlf_bottom.r_squared()
                check = False
            num += 1

        x = Symbol('x')

        # top
        for i in range(sec.pwlf_top.n_segments):
            eq = get_symbolic_eqn(sec.pwlf_top, i + 1)
            sec.top.append(eq)
            sec.f_top.append(lambdify([x], eq))
            sec.points_top = breaks1

        # bottom
        for i in range(sec.pwlf_bottom.n_segments):
            eq = get_symbolic_eqn(sec.pwlf_bottom, i + 1)
            sec.bottom.append(eq)
            sec.f_bottom.append(lambdify([x], eq))
            sec.points_bottom = breaks2
        secs.append(sec)
    return secs


def time_comp(secs):
    times = []
    for sec in secs:
        x_val = sec.points_top[-1]
        y_val = sec.f_top[-1](x_val)
        x, y = sym.symbols('x y')
        val = sec.t().subs({x: x_val, y: y_val})
        times.append(val)
    times2 = [times[0]]
    for time in times:
        add = True
        for i in range(len(times2)):
            if abs(time-times2[i]) < 0.01:
                add = False
        if add:
            times2.append(time)

    return times2


def pdf(vals, prob_times, dx, n_squares, index, j, test=False):
    times = list(np.arange(0, 10, 20*dx))
    a = list(sorted(vals))
    factor = 1/min(a)*min(prob_times)
    b = []
    for item in a:
        b.append(item*factor)
    cdf = [0]
    # compute cdf
    for t in times:
        num = cdf[-1]
        for i in range(num, len(b)):
            if b[i] <= t:
                num += 1
                continue
            else:
                cdf.append(num)
                break
    # compute pdf
    pdf = []
    for i in range(len(cdf) - 1):
        delta = (cdf[i+1] - cdf[i])/dx
        pdf.append(delta)
    fig, ax = plt.subplots(figsize=(10, 10))
    #print(f'length of inputs: {len(times)}, {len(pdf)}')
    ax.scatter(times, pdf, s=0.5)

    # plot discontinuities
    probs2 = []
    for item in prob_times:
        probs2.append(item)
        # probs2.append(factor*item)
    for t in probs2:
        if t > max(times):
            continue
        for i in range(len(times)):
            if t < times[i]:
                ax.scatter(t, pdf[i-1], s=20, color="red")
                break
    if test == True:
        print(prob_times)
    plt.show()
    plt.savefig(os.path.join(
        "results", f"{n_squares} - {index}", f"pdf - {j}"))
    plt.close(fig)
    return pdf


# Code from pwlf

x = Symbol('x')


def get_symbolic_eqn(pwlf_, segment_number):
    if pwlf_.degree < 1:
        raise ValueError('Degree must be at least 1')
    if segment_number < 1 or segment_number > pwlf_.n_segments:
        raise ValueError('segment_number not possible')
    # assemble degree = 1 first
    for line in range(segment_number):
        if line == 0:
            my_eqn = pwlf_.beta[0] + (pwlf_.beta[1])*(x-pwlf_.fit_breaks[0])
        else:
            my_eqn += (pwlf_.beta[line+1])*(x-pwlf_.fit_breaks[line])
    # assemble all other degrees
    if pwlf_.degree > 1:
        for k in range(2, pwlf_.degree + 1):
            for line in range(segment_number):
                beta_index = pwlf_.n_segments*(k-1) + line + 1
                my_eqn += (pwlf_.beta[beta_index]) * \
                    (x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()


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

# generate vectors for saddle connections on STS


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


def save_arrays_to_file(file_path, arrays_list):
    # Save arrays to a single NumPy file
    np.save(file_path, arrays_list)

# take files from saddle.py and load them into notebook


def load_arrays_from_file(file_path):
    # Load arrays from the NumPy file
    arrays_list = np.load(file_path, allow_pickle=True)

    # Ensure each element in the list is a NumPy array
    arrays_list = [np.array(array) for array in arrays_list]

    return arrays_list


def covolume(secs):
    sum = 0
    for i in range(len(secs)):

        all_points = []
        for item in secs[i].points_top:
            all_points.append(item)
        for item in secs[i].points_bottom:
            all_points.append(item)

        all_points = set(all_points)
        all_points = list(all_points)
        all_points.sort()

        m = 0
        n = 0
        for k in range(1, len(all_points)):
            if (not (all_points[k-1] >= secs[i].points_top[m] and all_points[k] <= secs[i].points_top[m+1])):
                m += 1
            if (not (all_points[k-1] >= secs[i].points_bottom[n] and all_points[k] <= secs[i].points_bottom[n+1])):
                n += 1
            top_eq = secs[i].f_top[m]
            bottom_eq = secs[i].f_bottom[n]
            upper = all_points[k]
            lower = all_points[k-1]
            x_ = float(secs[i].vec[0][0])
            y_ = float(secs[i].vec[1][0])

            # Perform the double integral
            def f(y, x):
                print(f"x: {x}, x_: {x_}, y: {y}, y_: {y_}")
                return y_ / (x * (x_ * x + y_ * y))
            result, error = integrate.dblquad(
                f, lower, upper, lambda x: bottom_eq(x), lambda x: top_eq(x))
            sum += result
    return sum


def read_df(n_squares, index, cusp):
    df = pd.read_csv(os.path.join(
        "results", f"{n_squares} - {index}", "df - " + str(cusp)))

    def read_vec(s):
        s = s.replace("[", "")
        s = s.replace("]", "")
        nums = s.split("\n")
        nums[0] = frac(nums[0])
        nums[1] = frac(nums[1])
        return np.array([[nums[0]], [nums[1]]])

    df["vec"] = df["vec"].apply(read_vec)
    return df


def winners1(vecs0, x_vals, m0, m1, y0, dx, dx_y):
    # dictionary for plotting
    saddle_dict = {}
    saddle_dict["x"] = []
    saddle_dict["y"] = []
    saddle_dict["lab"] = []
    saddle_dict["vec"] = []
    saddle_dict["time"] = []
    possible_vecs = []
    winners = []
    dz = 1 / 100
    vecs1 = np.hstack([arr.astype(float) for arr in vecs0])  # Ensure all are float64
    vecs = np.hstack(vecs0)  # Keep the original structure for output if needed

    # top edge
    # Stack list of numpy arrays into a single 2D array
    t0 = time()
    for a in np.arange(dz, 1, dz):
        # Matrix stays outside the inner loop
        Mab = np.array([[a, m0*a + 1/y0 - dx_y], [0, a]], dtype = 'float')

        new_vecs = Mab @ vecs1

        x_comps = new_vecs[0, :]
        y_comps = new_vecs[1, :]

        # Filter based on conditions (x > 0, y/x > 0, and x <= 1)
        valid_mask = (x_comps > 0) & (x_comps <= 1) & (y_comps / np.where(x_comps == 0,np.inf, x_comps) > 0)

        valid_x = x_comps[valid_mask]
        valid_y = y_comps[valid_mask]
        # Apply the mask to filter valid vectors
        valid_vecs = vecs[:, valid_mask]

        if valid_x.size == 0:
            winners.append(None)
            continue

        # Calculate slopes
        slopes = valid_y / valid_x

        # Find the minimum slope and handle continuity cases
        min_slope_idx = np.argmin(slopes)
        winner_slope = slopes[min_slope_idx]
        winner = valid_vecs[:, min_slope_idx]

        # Handle continuity by finding the smallest vector if slopes are close
        for i, slope in enumerate(slopes):
            if np.abs(slope - winner_slope) <= dx / 1000:
                if valid_vecs[0, i] < winner[0] or valid_vecs[1, i] < winner[1]:
                    winner = valid_vecs[:, i]

        winners.append(winner.reshape(2, 1))
    t1 = time()
    print("top done: " + str(t1 - t0))

    # verticals
    t0 = time()
    x_vals_side = np.arange(dx, 1, dx*100)
    for a in x_vals_side:
        y_vals = np.arange(m1*a + 1/y0 + dx_y, m0*a + 1/y0 - dx_y, dx_y*2)
        for b in y_vals:
            Mab = np.array([[a, b], [0, 1/a]], dtype = 'float')
            # Apply the transformation to all vectors at once
            new_vecs = Mab @ vecs1
    
            x_comps = new_vecs[0, :]
            y_comps = new_vecs[1, :]
    
            # Filter based on conditions (x > 0, y/x > 0, and x <= 1)
            valid_mask = (x_comps > 0) & (x_comps <= 1) & (y_comps / np.where(x_comps == 0,np.inf, x_comps) > 0)
    
            valid_x = x_comps[valid_mask]
            valid_y = y_comps[valid_mask]
            # Apply the mask to filter valid vectors
            valid_vecs = vecs[:, valid_mask]
    
            if valid_x.size == 0:
                winners.append(None)
                continue
    
            # Calculate slopes
            slopes = valid_y / valid_x
    
            # Find the minimum slope and handle continuity cases
            min_slope_idx = np.argmin(slopes)
            winner_slope = slopes[min_slope_idx]
            winner = valid_vecs[:, min_slope_idx]
    
            # Handle continuity by finding the smallest vector if slopes are close
            for i, slope in enumerate(slopes):
                if np.abs(slope - winner_slope) <= dx / 1000:
                    if valid_vecs[0, i] < winner[0] or valid_vecs[1, i] < winner[1]:
                        winner = valid_vecs[:, i]
            winners.append(winner.reshape(2, 1))
    t1 = time()
    print("verticals done: " + str(t1 - t0))

    # side edge
    t0 = time()
    y_vals = np.arange(m1*a + 1/y0 + dx, m0*a + 1/y0 - dx, dx)
    for b in y_vals:
        Mab = np.array([[1-dx, b], [0, 1/(1-dx)]], dtype = 'float')
        # Apply the transformation to all vectors at once
        new_vecs = Mab @ vecs1

        x_comps = new_vecs[0, :]
        y_comps = new_vecs[1, :]

        # Filter based on conditions (x > 0, y/x > 0, and x <= 1)
        valid_mask = (x_comps > 0) & (x_comps <= 1) & (y_comps / np.where(x_comps == 0,np.inf, x_comps) > 0)

        valid_x = x_comps[valid_mask]
        valid_y = y_comps[valid_mask]
        # Apply the mask to filter valid vectors
        valid_vecs = vecs[:, valid_mask]

        if valid_x.size == 0:
            winners.append(None)
            continue

        # Calculate slopes
        slopes = valid_y / valid_x

        # Find the minimum slope and handle continuity cases
        min_slope_idx = np.argmin(slopes)
        winner_slope = slopes[min_slope_idx]
        winner = valid_vecs[:, min_slope_idx]

        # Handle continuity by finding the smallest vector if slopes are close
        for i, slope in enumerate(slopes):
            if np.abs(slope - winner_slope) <= dx / 1000:
                if valid_vecs[0, i] < winner[0] or valid_vecs[1, i] < winner[1]:
                    winner = valid_vecs[:, i]
        winners.append(winner.reshape(2, 1))
    t1 = time()
    print("side done: " + str(t1 - t0))

    # diagonal
    t0 = time()
    for a in np.arange(0 + dz, 1, dz):
        Mab = np.array([[a, m1*a + 1/y0 + dx_y], [0, a]], dtype = 'float')
        # Apply the transformation to all vectors at once
        new_vecs = Mab @ vecs1

        x_comps = new_vecs[0, :]
        y_comps = new_vecs[1, :]

        # Filter based on conditions (x > 0, y/x > 0, and x <= 1)
        valid_mask = (x_comps > 0) & (x_comps <= 1) & (y_comps / np.where(x_comps == 0,np.inf, x_comps) > 0)

        valid_x = x_comps[valid_mask]
        valid_y = y_comps[valid_mask]
        # Apply the mask to filter valid vectors
        valid_vecs = vecs[:, valid_mask]

        if valid_x.size == 0:
            winners.append(None)
            continue

        # Calculate slopes
        slopes = valid_y / valid_x

        # Find the minimum slope and handle continuity cases
        min_slope_idx = np.argmin(slopes)
        winner_slope = slopes[min_slope_idx]
        winner = valid_vecs[:, min_slope_idx]

        # Handle continuity by finding the smallest vector if slopes are close
        for i, slope in enumerate(slopes):
            if np.abs(slope - winner_slope) <= dx / 1000:
                if valid_vecs[0, i] < winner[0] or valid_vecs[1, i] < winner[1]:
                    winner = valid_vecs[:, i]
        winners.append(winner.reshape(2, 1))
    t1 = time()
    print("diagonal done: " + str(t1 - t0))
    
    winners2 = []
    for winner in winners:
        try:
            if winner == None:
                continue
        except:
            winners2.append(winner)
    # Flatten winners and store as a list of lists
    flattened_winners = [winner.flatten().tolist() for winner in winners2]
    
    # Find unique entries
    unique_winners = []
    for item in flattened_winners:
        if item not in unique_winners:
            unique_winners.append(item)
    
    # Manually assemble each unique winner into a 2D column vector
    possible_vecs = []
    for unique_item in unique_winners:
        # Manually create the column vector
        vector = np.array([[unique_item[0]], [unique_item[1]]], dtype=object)
        possible_vecs.append(vector)

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

    #print(possible_vecs)
    #print(possible_vecs_float)

    for i in range(len(possible_vecs_float)):
        # apply Mab matrix, perform horocycle flow and find time t to horizontal
        a = horo@(Mab@possible_vecs_float[i])
        #print(a)
        t_dict[i] = lambdify([x, y], solve(a[1][0], t)[0])

    # for each point (a,b) in the poincare section, apply the Mab matrix to each vector and look for "winners". Winners have smallest possible slope that is greater than zero and 0 < x-component <= 1
    for a in x_vals:
        y_vals = np.arange(m1*a + 1/y0 + dx_y, m0*a + 1/y0 - dx_y, dx_y)
        for b in y_vals:
            check = 0
            winner_slope = None
            winner = None
            Mab = np.array([[a, b], [0, 1/a]], dtype = 'float')
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
                saddle_dict["lab"].append(len(vecs0))
                saddle_dict["time"].append(1000)

    df = pd.DataFrame.from_dict(saddle_dict)
    return df


def sec_setup2(df, dx_y):
    sec_list = []
    global labs
    labs = df["lab"].unique()
    vec_order = []
    vec_dict = {}
    #print(labs)
    for lab in labs:
        sec_dict = {}
        df1 = df[df["lab"] == lab]
        vec_order.append(df1['vec'].iloc[int(0)])
        vec_dict[lab] = df1['vec'].iloc[int(0)]
        xs = df1["x"]
        xs = sorted(list(set(xs.tolist())))
        y_tops = []
        y_bottoms = []
        for x in xs:
            # for a given "x" find the max and minimum y-values
            y_top = max(df1[df1["x"] == x]["y"])
            y_bottom = min(df1[df1["x"] == x]["y"])
            # ensures the section is convex, not concave
            if len(df1[df1["x"] == x]["y"]) < (y_top - y_bottom)/dx_y:
                print("len: " + str(len(df1[df1["x"] == x]["y"])))
                print("ytop: " + str(y_top))
                print("ybottom: " + str(y_bottom))
                print("dx_y: " + str(dx_y))
                print(x)
                print(df1[df1["x"] == x]["y"])
                raise ValueError(
                    "Section has more than 2 points for a given 'x'")
            y_tops.append(y_top)
            y_bottoms.append(y_bottom)
        y_tops = np.array(y_tops)
        y_bottoms = np.array(y_bottoms)
        xs = np.array(xs)
        sec_dict['x'] = xs
        sec_dict['top'] = y_tops
        sec_dict['bottom'] = y_bottoms
        sec_list.append(sec_dict)
    return sec_list, vec_order, vec_dict


def sec_comp2(df, sec_list, vec_order, vec_dict, dx, dx_y, m1, y0):
    from fractions import Fraction as frac
    x = Symbol('x')
    
    secs = []
    problem_xs = []
    for i in range(len(sec_list)):
        xs = sec_list[i]['x']
        problem_xs.append(xs[0])
    problem_xs = sorted(list(set(problem_xs)))
    #print(problem_xs)
    problems = [problem_xs[0]]
    for i in range(len(problem_xs) - 1):
        if abs(problem_xs[i] - problem_xs[i+1]) <= 0.01:
            continue
        problems.append(problem_xs[i+1])
    problems.append(1)
    #print(problems)
    prob_mids = []
    for i in range(len(problems) - 1):
        prob_mids.append((problems[i] + problems[i+1])/2)
    i = 0
    use_points = []
    for item in np.arange(dx, 1, dx):
        if item < prob_mids[i]:
            continue
        i += 1
        use_points.append(item)
        if i == len(prob_mids):
            break
    #print(use_points)
        
    for i in range(len(sec_list)):
        xs = sec_list[i]['x']
        top = sec_list[i]['top']
        bottom = sec_list[i]['bottom']
        sec = Section(x, top, bottom)
        sec.vec = vec_order[int(i)]

        bottom_lab_list = []
        for x_, y_ in zip(xs, bottom):
            if x_ not in use_points:
                continue
            poss_ys = sorted(list(df[(df["x"] == x_)]["y"]), reverse = True)
            y_index = poss_ys.index(y_)
        
            # Check if index + 1 is valid
            if y_index + 1 < len(poss_ys):
                # Get the value at index + 1
                next_y = poss_ys[y_index + 1]
                filtered_df = df[(df["x"] == x_) & (df["y"] == next_y)]
                lab_value = filtered_df["lab"].iloc[0]  # Get the value in the "lab" column of the first row
                if lab_value not in bottom_lab_list:
                    bottom_lab_list.append(lab_value)
            else:
                if -1 not in bottom_lab_list:
                    bottom_lab_list.append(int(-1))

        eqs = []
        #print(bottom_lab_list)
        for lab in bottom_lab_list:
            if lab != -1:
                vec = vec_dict[lab]
                eqs.append(-frac(vec[0][0],vec[1][0]) * x + frac(int(1)/vec[1][0]))
            else:
                eqs.append(frac(m1) * x + frac(int(1)/y0))

        points = [xs[0]]
        for i in range(len(eqs) - 1):
            eq1 = eqs[i]
            eq2 = eqs[i + 1]
            # Set eq1 = eq2 and solve for x
            intersection_x = solve(Eq(eq1, eq2), x)
            if intersection_x:  # If there is a solution
                points.append(intersection_x[0])
            else:
                print(f"intersection between bottom equation {i} and bottom equation {i+1} point failed")

        sec.top.append(-frac(sec.vec[0][0], sec.vec[1][0]) * x + frac(int(1)/sec.vec[1][0]))
        sec.f_top.append(lambdify([x], sec.top[0]))

        # modifying the first entry in points so that it is not the rounded version. Finding the x-componennt of the intersection between the    
        # "top" equation and the first "bottom" equation
        eq1 = sec.top[0]
        eq2 = eqs[0]

        intersection_x = solve(Eq(eq1, eq2), x)
            if intersection_x:  # If there is a solution
                points[0] = intersection_x[0]
            else:
                print("intersection between top and first bottom point failed")
        
        sec.points_top = [points[0], points[-1]]
        
        for eq in eqs:
            sec.bottom.append(eq)
            sec.f_bottom.append(lambdify([x], eq))
            sec.points_bottom = points
            
        secs.append(sec)
    return secs

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
