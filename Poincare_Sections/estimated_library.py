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

class Section_estimated:
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

    # output time equation
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

# Code from pwlf, used to get piecewise representation of the sub-sections
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

# this is only used for computing the estimated pdfs for a given poincare section


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

# this is only used for computing the estimated pdfs for a given poincare section


def sec_comp(sec_list, dx):
    secs = []
    for i in range(len(sec_list)):
        x = sec_list[i]['x']
        top = sec_list[i]['top']
        bottom = sec_list[i]['bottom']
        sec = Section_estimated(x, top, bottom)
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
    # print(f'length of inputs: {len(times)}, {len(pdf)}')
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
        "results", f"{n_squares}_{index}", f"pdf_{j}"))
    plt.close(fig)
    return pdf