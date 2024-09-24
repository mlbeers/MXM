import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import numpy as np
from matplotlib import pyplot as plt
import os
import pwlf
import sympy as sym
from sympy import Symbol
from sympy import solve, lambdify
import os
import math
from surface_dynamics.all import *
from time import time
from scipy import integrate

def generators(perm, vecs0):
    #find the generators of each cusp of the STS
    generators = []
    a = perm.veech_group().cusps()
    for item in a:
        m = perm.veech_group().cusp_data(item)[0]
        generators.append(m.matrix())
    return generators
    
def poincare_details(perm, vecs0, generators):
    #find the generators of each cusp of the STS
    generators = []
    a = perm.veech_group().cusps()
    for item in a:
        m = perm.veech_group().cusp_data(item)[0]
        generators.append(m.matrix())
    #find the eigenvectors for each generator and make sure they are 1
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
    #find the eigenvectors for each generator
    eigenvecs = []
    for matrix in generators:
        vec = matrix.eigenvectors_right()[0][1][0]
        vec = np.array([[vec[0]],[vec[1]]])
        eigenvecs.append(vec)
    #find the magnitude, slope, x-direction, and y-direction of each eigenvector
    saddle_vecs = []
    for vec in eigenvecs:
        # Update: Try to optimize the original 44-93 lines
        mag_vec = np.linalg.norm(vec) #  use numpy's linalg.norm to calc vector's norm
        slope_vec = float('inf') if vec[0] == 0 else vec[1] / vec[0]
        x_sign_vec = np.sign(vec[0])
        y_sign_vec = np.sign(vec[1])

        saddle_vec = None
        min_mag = float('inf')

        for saddle in vecs0:
            mag_saddle = np.linalg.norm(saddle)
            slope_saddle = float('inf') if saddle[0] == 0 else saddle[1] / saddle[0]
            x_sign_saddle = np.sign(saddle[0])
            y_sign_saddle = np.sign(saddle[1])

            # Check if the slope and direction match
            if (slope_vec == slope_saddle and x_sign_vec == x_sign_saddle and y_sign_vec == y_sign_saddle):

                # Update the smallest saddle connection
                if mag_saddle < min_mag:
                    min_mag = mag_saddle
                    saddle_vec = saddle

        if saddle_vec is None:
            raise ValueError(f"No saddle vec for eigenvector {vec}")
        saddle_vecs.append(saddle_vec)

    #find the counter-clockwise angle from the x-axis to the eigenvectors
    thetas = []
    for i in range(len(saddle_vecs)):
        mag = (saddle_vecs[i][0]**2 + saddle_vecs[i][1]**2)**0.5
        theta = np.arccos(np.dot(np.array([[1,0]]),saddle_vecs[i])/mag)
        if saddle_vecs[i][1] < 0:
            theta = 2 * math.pi - theta
        thetas.append(theta)
    #find the rotation matrix that takes the vector (1,0) to the vector in the direction of each eigenvector
    rots = []
    for theta in thetas:
        rot = np.array([[round(np.cos(theta)[0][0],5),round(-np.sin(theta)[0][0],5)],[round(np.sin(theta)[0][0],5), round(np.cos(theta)[0][0],5)]])
        rots.append(rot)

    #find a constant value such that mult*rot@(1,0) = saddle_vec while accounting for rounding errors and zero matrix inputs
    mults = []
    for i in range(len(rots)):
        matrix = rots[i]@np.array([[1],[0]])
        if matrix[0][0] != 0:
            mult1 = saddle_vecs[i][0][0]/matrix[0][0]
        else:
            mult1 = 0
        if matrix[1][0] != 0:
            mult2 = saddle_vecs[i][1][0]/matrix[1][0]
        else:
            mult2 = 0
        if mult1 != 0 and mult2 == 0:
            mult = mult1
        elif mult1 == 0 and mult2 != 0:
            mult = mult2
        elif mult1 == 0 and mult2 == 0:
            raise ValueError('both mults equal zero')
        elif abs(mult1 - mult2) <= 0.001:
            mult = mult1
        elif abs(mult1 - mult2) >= 0.001:
            raise ValueError(f'mults are different {mult1}, {mult2}')
        mults.append(mult)
        mult1 = None
        mult2 = None
        mult = None
    #find c_inv and c
    Cs = []
    C_invs = []
    for i in range(len(mults)):
        c_inv = mults[i]*rots[i]
        c = np.linalg.inv(c_inv)
        Cs.append(c)
        C_invs.append(c_inv)
    #alpha is the top right value of the matrix M = c @ generator @ c_inv. M must have 1s on the diagonal and 0 in the bottom left
    alphas = []
    Ms = []
    for i in range(len(generators)):
        M = Cs[i]@generators[i]@C_invs[i]
        Ms.append(M)
        if M[1][0] >= 1/1000000 and M[1][0] <= -1/1000000:
            raise ValueError(f"Wrong conjugate matrix\nC: {Cs[i]}\nC_inv: {C_invs[i]}\nM: {M}\ngenerator: {generators[i]}")
        alphas.append(round(M[0][1], 5))
    return alphas, Cs, C_invs, eigs, Ms, generators, eigenvecs
    
def setup(alpha, c, eig, vecs0, dx, improved = True):
    x_vals = np.arange(dx, 1, dx)
    #for the poincare section with matrix c, the original saddle connections are acted on by this matrix
    vecs1 = c@vecs0
    #find vectors such that -10 <= x <= 10 and 0 < y <= 10
    vecs = []
    for item in vecs1:
        vecs.append(item)

    #section vector is vector with smallest y-component with y > 0 and corresponding x-component with 0 <= x <= 1 with smallest slope. It defines the y-intercept of the section and partially defines slopes of the lines of the section
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

    #slopes of top and bottom lines of section
    if improved == True:
        m0 = -x0/y0
        m1 = -(x0/y0 + alpha)
    else:
        m0 = 0
        m1 = -alpha
        x0 = 0
        y0 = 1

    #defines vertical step when finding winners
    global dx_y
    dx_y = (m0 - m1) * dx
    
    return vecs, x_vals, m0, m1, x0, y0, dx_y, z


def winners(vecs, x_vals, m0, m1, y0, dx, dx_y):
    #dictionary for plotting
    saddle_dict = {}
    saddle_dict["x"] = []
    saddle_dict["y"] = []
    saddle_dict["lab"] = []
    saddle_dict["vec"] = []
    saddle_dict["time"] = []
    possible_vecs = []
    winners = []

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
        #if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                elif abs(y/x - winner_slope) <= dx/1000:
                    if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                        winner = vec
                        continue
                elif y/x < winner_slope:
                    winner_slope = y/x
                    winner = vec
                    continue
        winners.append(winner)

    #diagonal
    for a in np.arange(0+ dz, 1, dz):
        check = 0
        winner_slope = None
        winner = None
        Mab = np.array([[a, m1*a + 1-dx/y0 + dx_y], [0, a]])
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
       #if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                elif abs(y/x - winner_slope) <= dx/1000:
                    if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                        winner = vec
                        continue
                elif y/x < winner_slope:
                    winner_slope = y/x
                    winner = vec
                    continue
        winners.append(winner)

    #side edge
    y_vals = np.arange(m1 + (1-dx)/y0 + dx_y, m0 + (1-dx)/y0 - dx_y, dz*(m0-m1))
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
        #if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
                elif abs(y/x - winner_slope) <= dx/1000:
                    if vec[0][0] < winner[0][0] or vec[1][0] < winner[1][0]:
                        winner = vec
                        continue
                elif y/x < winner_slope:
                    winner_slope = y/x
                    winner = vec
                    continue
        winners.append(winner)

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
    
    #dictionary for vector labels
    for i in range(len(possible_vecs)):
        label_dict[i] = possible_vecs[i]

    #for each vector, there is a time function defined as f(a,b) where a,b are points in the poincare section
    global t_dict
    t_dict = {}

    x, y, t = sym.symbols('x y t')
    Mab = np.array([[x, y], [0, 1/x]])
    horo = np.array([[1, 0], [-t, 1]])

    for i in range(len(possible_vecs)):
        #apply Mab matrix, perform horocycle flow and find time t to horizontal
        a = horo@(Mab@possible_vecs[i])
        t_dict[i] = lambdify([x,y], solve(a[1][0], t)[0])
    
    #for each point (a,b) in the poincare section, apply the Mab matrix to each vector and look for "winners". Winners have smallest possible slope that is greater than zero and 0 < x-component <= 1
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
            #if you have two potential winners like (m,n) and 2*(m,n), make (m,n) winner for continuity and plotting purposes
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
                    saddle_dict["time"].append(t_dict[i](a,b))
                    if saddle_dict["time"][-1] < 0:
                        saddle_dict["time"][-1] = 1000
            #if there is no winner, at (a,b), add a label so the df and plot can still be made. These section will later be made blank for trouble-shoooting
            if check == 0:
                saddle_dict["lab"].append(len(vecs))
                saddle_dict["time"].append(1000)

    df = pd.DataFrame.from_dict(saddle_dict)
    return df

def plot(df, vecs, c, j, n_squares, index, test = False):
    fig, ax = plt.subplots(figsize=(10, 10))
    #plot winners
    ax.scatter(df[df["lab"] != len(vecs)]["x"], df[df["lab"] != len(vecs)]["y"], c = df[df["lab"] != len(vecs)]["lab"], s = 0.1)
    #for points with no winner, make white
    ax.scatter(df[df["lab"] == len(vecs)]["x"], df[df["lab"] == len(vecs)]["y"], c = "white", alpha = 0.5, s = 0.1)
    ax.set_aspect('auto', adjustable='box')

    #plot outline of poincare section
    min_x = min(df["x"])
    y1 = max(df[df["x"] == min_x]["y"])
    y2 = min(df[df["x"] == min_x]["y"])
    max_x = max(df["x"])
    y3 = max(df[df["x"] == max_x]["y"])
    y4 = min(df[df["x"] == max_x]["y"])
    ax.plot([min_x, max_x, max_x, min_x, min_x], [y1, y3, y4, y2, y1], c = "black")
    
    plt.title(str(c) + "\n", loc = "left", fontsize = 25)
    if test == False:
        plt.savefig(os.path.join("results", f"{n_squares} - {index}", "section - " + str(j)))
    #for troubleshooting
    if test == True:
        #display vectors on the right edge of the section from top to bottom
        labs = list(df[df["x"] == max(df["x"])]["lab"].unique())
        labs.reverse()
        output = []
        for lab in labs:
            if lab == len(vecs):
                continue
            output.append([label_dict[lab][0][0], label_dict[lab][1][0]])
        print(output)
    if len(df[df["lab"] == len(vecs)]) != 0:
        raise ValueError("Poincare section has empty portion")
    plt.show()
    plt.close(fig)
    
class Section:
    def __init__(self, x, top, bottom):
        self.vec = None
        self.pwlf_top = pwlf.PiecewiseLinFit(x, top)
        self.pwlf_bottom = pwlf.PiecewiseLinFit(x, bottom)
        #equations
        self.top = []
        self.bottom = []
        #lambdified equations
        self.f_top = []
        self.f_bottom = []
        #points
        self.points_top = None
        self.points_bottom = None
        
    def t(self):
        x, y, t = sym.symbols('x y t')
        Mab = np.array([[x, y], [0, 1/x]])
        horo = np.array([[1, 0], [-t, 1]])
        a = horo@(Mab@self.vec)
        return solve(a[1][0], t)[0]

    #find the time equation in terms of y
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
            #for a given "x" find the max and minimum y-values
            y_top = max(df1[df1["x"] == x]["y"])
            y_bottom = min(df1[df1["x"] == x]["y"])
            #ensures the section is convex, not concave
            if len(df1[df1["x"] == x]["y"]) < (y_top - y_bottom)/dx_y:
                print("len: " + str(len(df1[df1["x"] == x]["y"])))
                print("ytop: " + str(y_top))
                print("ybottom: " + str(y_bottom))
                print("dx_y: " + str(dx_y))
                print(x)
                print(df1[df1["x"] == x]["y"])
                raise ValueError("Section has more than 2 points for a given 'x'")
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

        #use piece-wise linear regression to find the equations of the lines for subsection
        #top
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

        #bottom
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

        #top
        for i in range(sec.pwlf_top.n_segments):
            eq = get_symbolic_eqn(sec.pwlf_top, i + 1)
            sec.top.append(eq)
            sec.f_top.append(lambdify([x], eq))
            sec.points_top = breaks1

        #bottom
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

def pdf(vals, prob_times, dx, n_squares, index, j, test = False):
    times = list(np.arange(0, 10, 20*dx))
    a = list(sorted(vals))
    factor = 1/min(a)*min(prob_times)
    b = []
    for item in a:
        b.append(item*factor)
    cdf = [0]
    #compute cdf
    for t in times:
        num = cdf[-1]
        for i in range(num, len(b)):
            if b[i] <= t:
                num += 1
                continue
            else:
                cdf.append(num)
                break
    #compute pdf
    pdf = []
    for i in range(len(cdf) - 1):
        delta = (cdf[i+1] - cdf[i])/dx
        pdf.append(delta)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(times, pdf, s = 0.5)
        
    #plot discontinuities
    probs2 = []
    for item in prob_times:
        probs2.append(item)
        #probs2.append(factor*item)
    for t in probs2:
        if t > max(times):
            continue
        for i in range(len(times)):
            if t < times[i]:
                ax.scatter(t, pdf[i-1], s = 20, color = "red")
                break
    if test == True:
        print(prob_times)
    plt.show()
    plt.savefig(os.path.join("results", f"{n_squares} - {index}", f"pdf - {j}"))
    plt.close(fig)
    return pdf


#Code from pwlf

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
                my_eqn += (pwlf_.beta[beta_index])*(x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()

#Code from Sunrose 
#Gives non-visibility tori for testing
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
            count_obstructed+= item.teichmueller_curve().orbit_graph().num_verts()
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
def vectors2(perm, length = 200):
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
        item = np.array([[vec[0]],[vec[1]]])
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
        for k in range (1, len(all_points)):
            if (not(all_points[k-1] >= secs[i].points_top[m] and all_points[k] <= secs[i].points_top[m+1])):
                m += 1
            if (not(all_points[k-1] >= secs[i].points_bottom[n] and all_points[k] <= secs[i].points_bottom[n+1])):
                n += 1
            top_eq = secs[i].f_top[m]
            bottom_eq = secs[i].f_bottom[n]
            upper = all_points[k]
            lower = all_points[k-1]
            
            # Perform the double integral
            f = lambda y,x :secs[i].vec[1][0] / (x * (secs[i].vec[0][0] * x + secs[i].vec[1][0] * y))
            result, error = integrate.dblquad(f, lower, upper, lambda x: bottom_eq(x), lambda x: top_eq(x)) 
            sum += result
    return sum

def read_df(n_squares, index, cusp):
    df = pd.read_csv(os.path.join("results", f"{n_squares} - {index}", "df - " + str(cusp)))

    def read_vec(s):
        s = s.replace("[", "")
        s = s.replace("]", "")
        nums = s.split("\n")
        nums[0] = float(nums[0])
        nums[1] = float(nums[1])
        return np.array([[nums[0]],[nums[1]]])

    df["vec"] = df["vec"].apply(read_vec)
    return df
