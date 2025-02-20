import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from flatsurf import *
import os
import pwlf
from surface_dynamics.all import *
from Library import *
from Library import Section
import math
from time import time
import copy
from scipy import integrate
import sympy as sym
from sympy import Symbol, solve, lambdify
import traceback
import dill
import sys
import unittest
from surface_dynamics.all import Origami
from utils import load_arrays_from_file  # testing
import re
import subprocess
from sage.all import *
import numpy as np
from fractions import Fraction as frac
M = mathematica
from IPython.display import display, Math
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy import symbols, And, Piecewise, Add, N
from sympy.core.relational import Relational
from integration_functions import *

t = sp.Symbol('t')

def run_integrals(n_squares, index, a):
    boundary_points = set()
    list_functions = []
    M('ClearAll["Global`*"]')
    allFunctions = M('allFunctions= {}')
    totalVol = M('totalVol = 0')
    
    for j in range(len(a[0])):
        df = read_df(n_squares, index, j)
        print(f"section {j}")
    
        with open(os.path.join("results", f"{n_squares} - {index}", f"secs_integrals - {j}.dill"), 'rb') as f:
            secs2 = dill.load(f)  
    
        sections = get_secs(secs2)
    
        # Create directory for storing LaTeX outputs
        results_dir = os.path.join("results", f"{n_squares} - {index}", f"cusp - {j}")
        os.makedirs(results_dir, exist_ok=True)
    
        # LaTeX file path
        latex_file_path = os.path.join(results_dir, f"equations_{j}.tex")
        pdf_file_path = os.path.join(results_dir, f"equations_{j}.pdf")
    
        with open(latex_file_path, "w") as latex_file:
            # Start LaTeX document
            latex_file.write("\\documentclass{article}\n")
            latex_file.write("\\usepackage{amsmath}\n")
            latex_file.write("\\begin{document}\n\n")
    
            for i, section in enumerate(sections):
                print(section)
                if len(section) == 5:
                    x0, y0, top, bottom1, left = section
                    bottoms = [bottom1]
                    points = []
                elif len(section) == 7:
                    x0, y0, top, bottom1, bottom2, point1, left = section
                    bottoms = [bottom1, bottom2]
                    points = [point1]
                elif len(section) == 9:
                    x0, y0, top, bottom1, bottom2, bottom3, point1, point2, left = section
                    bottoms = [bottom1, bottom2, bottom3]
                    points = [point1, point2]
    
                eqs, totalVol = integrals(x0, y0, top, bottoms, points, left)
    
                # Convert equation to LaTeX format
                expression_string = repr(eqs[1])
                es, cs = parse_piecewise_sp(expression_string)
                for c in cs:
                    for c_ in c:
                        boundary_points.add(c_)
    
                eq_tuples = []
                for eq, cond in zip(es, cs):
                    eq_tuples.append((
                        parse_expr(eq, local_dict={'t': t, 'sp': sp}),
                        sp.And(sp.simplify(cond[0]) <= t, t < sp.simplify(cond[1]))
                    ))
    
                f = sp.Piecewise(*eq_tuples)
                list_functions.append(f)
                latex_expr = sp.latex(f)
                latex_expr = latex_expr.replace(r"\text{for}\: t", "").replace(r"\geq", "").replace(r"\wedge", r"\leq")
    
                # Write equation to LaTeX file
                latex_file.write(f"Equation {i}:\n\\[\n{latex_expr}\n\\]\n\n")
    
                # Generate graph for visualization
                graph_piece(f, cs, n_squares, index, j, i)
    
            # End LaTeX document
            latex_file.write("\n\\end{document}\n")
    
        print(f"LaTeX file saved: {j}")
    
        # Compile the LaTeX file into a PDF while suppressing output
        try:
            with open(os.devnull, "w") as FNULL:
                subprocess.run(
                    ["pdflatex", "-output-directory", results_dir, latex_file_path],
                    stdout=FNULL, stderr=FNULL, check=True
                )
    
            print(f"PDF successfully created: {j}")
    
            # Remove auxiliary files (.log, .aux)
            for ext in [".log", ".aux"]:
                aux_file = os.path.join(results_dir, f"equations_{j}{ext}")
                if os.path.exists(aux_file):
                    os.remove(aux_file)
    
        except subprocess.CalledProcessError:
            print(f"Error: Failed to compile {j} into a PDF.")

    # Open the file in write mode
    totalVol = M('totalVol = FullSimplify[totalVol]')
    totalVol2 = M('totalVolN = N[totalVol]')
    target = M('target = N[(Pi^2/6)*54]')
    save_path = os.path.join("results", f"{n_squares} - {index}", "coVolume.txt")
    with open(save_path, "w") as file:
        file.write(f"{totalVol}\n\n\n")
        file.write(f"rounded coVol: {totalVol2}\n")
        file.write(f"rounded target: {target}\n")
    boundary_points = sorted([sp.simplify(bp) for bp in boundary_points])
    return list_functions, boundary_points

def create_combined_piecewise(piecewise_list, boundary_points):
    combined_piecewise = []

    for i in range(len(boundary_points) - 1):
        lower_bound = boundary_points[i]
        upper_bound = boundary_points[i + 1]

        # Create the condition for the current interval
        interval_condition = (t >= lower_bound) & (t < upper_bound)

        # Find all equations active in this interval
        active_equations = []
        for pw in piecewise_list:
            for expr, cond in pw.args:
                if isinstance(cond, And):
                    # Check if the interval is fully contained within the condition
                    if sp.simplify(And(interval_condition, cond)) == interval_condition:
                        active_equations.append(expr)
                elif isinstance(cond, Relational):
                    print("yup")
                    # Check if the interval is fully contained within the condition
                    if sp.simplify(And(interval_condition, cond)) == interval_condition:
                        active_equations.append(expr)

        # Add the active equations together
        if active_equations:
            combined_expression = Add(*active_equations)
            combined_piecewise.append((combined_expression, interval_condition))

    # Create the final combined Piecewise function
    return Piecewise(*combined_piecewise)

def get_secs(secs):
    outputs = []
    for sec in secs:
        output = []
        output.append("x0 = " + str(sec.vec[0][0]))
        output.append("y0 = " + str(sec.vec[1][0]))
        output.append("top = " + str(sec.top[0]))
        for i in range(len(sec.bottom)):
            output.append("bottom" + str(i+1) + " = " + str(sec.bottom[i]))
        for i in range(len(sec.points_bottom)):
            if i == 0:
                continue
            output.append("point" + str(i) + " = " + str(sec.points_bottom[i]))
        output.append("left = " + str(simplify(sec.points_bottom[0], 20)))
        outputs.append(output)
    return outputs

def simplify_eq(expr, den):
    dict = expr.as_coefficients_dict()
    b = dict[1]
    m = dict[x]

    m1 = None
    b1 = None

    m_num = None
    m_den = None
    b_num = None
    b_den = None
    
    m_diff = 1
    b_diff = 1
    m_decimal = m - int(m)
    b_decimal = b - int(b)
    
    for i in range(den, 0, -1):
        for j in range(i+1):
            if(m < 0):
                m_result = abs(j/i + m_decimal)
            else:
                m_result = abs(j/i - m_decimal)

            if(b < 0):
                b_result = abs(j/i + b_decimal)
            else:
                b_result = abs(j/i - b_decimal)

            if(m_result <= m_diff):
                m_diff = m_result
                if (m < 0):
                    m_num = int(m)*i - j
                    m_den = i
                else:
                    m_num = int(m)*i + j
                    m_den = i

            if(b_result <= b_diff):
                b_diff = b_result
                if (b < 0):
                    b_num = int(b)*i - j
                    b_den = i
                else:
                    b_num = int(b)*i + j
                    b_den = i
                    
    if(m_den == 1):
        m1 = str(m_num)
    else:
        m1 = str(m_num) + "/" + str(m_den)

    if(b_den == 1):
        b1 = str(b_num)
    else:
        b1 = str(b_num) + "/" + str(b_den)
        
    return  m1 + "*x + " + b1

def simplify(m, den):

    m1 = None

    m_num = None
    m_den = None
    
    m_diff = 1
    m_decimal = m - int(m)
    
    for i in range(den, 0, -1):
        for j in range(i+1):
            if(m < 0):
                m_result = abs(j/i + m_decimal)
            else:
                m_result = abs(j/i - m_decimal)

            if(m_result <= m_diff):
                m_diff = m_result
                if (m < 0):
                    m_num = int(m)*i - j
                    m_den = i
                else:
                    m_num = int(m)*i + j
                    m_den = i
                    
    if(m_den == 1):
        m1 = str(m_num)
    else:
        m1 = str(m_num) + "/" + str(m_den)
        
    return  m1

def nb_1(x0, y0, top, bottom1, left):
    M('ClearAll["Global`*"]')
    x0 = M(f'x0 = {x0}')
    y0 = M(f'y0 = {y0}')
    top = M(f'top = {top}')
    bottom1 = M(f'bottom1 = {bottom1}')
    left = M(f'left = {left}')
    func = M('func = 1/(t*x) - x0/y0*x')
    timeEnter = M('timeEnter = t /. Solve[y0 - t == 0, t][[1]]')
    symbolicZeros = M('symbolicZeros = Solve[bottom1 - func == 0, x]')
    r1 = M('r1 = x /. symbolicZeros[[1]]')
    r2 = M('r2 = x /. symbolicZeros[[2]]')
    timeEdge = M(' timeEdge = t /. Solve[r2 - r1 == 0, t][[1]]')
    try:
        tR = M('tR = t /. Solve[r2 - 1 == 0, t][[1]]')
    except:
        tR = None
    
    try:
        tL = M('tL = t /. Solve[r1 - left == 0, t][[1]]')
    except:
        tL = None
    
    try:
        tRA = M('tRA = t /. Solve[r1 - 1 == 0, t][[1]]')
    except:
        tRA = None
    
    try:
        tLA = M('tLA = t /. Solve[r2 - left == 0, t][[1]]')
    except:
        tLA = None
    
    #print(f"time Right End = {tR}")
    #print(f"time Left End = {tL}")
    #print()
    #print(f"Right Alt = {tRA}")
    #print(f"Left Alt = {tLA}")
    
    if tR is None and tRA is None:
        tR = M('tR = 20')
    if tL is None and tLA is None:
        tL = M('tL = 20')
    
    if tR is None and tRA is not None:
        tR = M('tR = tRA')
        if timeEdge >= tR:
            timeEdge = M('timeEdge = tR')
    if tL is None and tLA is not None:
        tL = M('tL = tLA')
        if timeEdge >= tL:
            timeEdge = M('timeEdge = tL')

    #print()
    #print(f"new Right End = {tR}")
    #print(f"new Left End = {tL}")
    #print(f"Edge = {timeEdge}")

    try:
        f1 = M('f1 = Simplify[D[Integrate[1, {x, y0/t, 1}, {y, func, top}], t]]')
    except:
        f1 = -1
    try:
        f2 = M('f2 = Simplify[f1 - Simplify[D[Integrate[1, {x, r1, r2}, {y, func, bottom1}], t]]]')
    except: 
        f2 = -1
    try:
        f3 = M('f3 = Simplify[f1 - Simplify[D[Integrate[1, {x, r1, 1}, {y, func, bottom1}], t]]]')
    except:
        f3 = -1
    try:
        fL = M('fL = Simplify[D[Integrate[1, {x, left, 1}, {y, func, top}], t]]')
    except:
        fL = -1
    try:
        f4 = M('f4 = Simplify[fL - Simplify[D[Integrate[1, {x, left, r2}, {y, func, bottom1}], t]]]')
    except:
        f4 = -1
    
    combined = M('combined = Piecewise[{{0, 0 <= t < timeEnter},{f1, timeEnter <= t < timeEdge},{f2, timeEdge <= t < Min[tR, tL]},{f3, tR <= t < tL},{f4, tL <= t < tR},{0, Max[tR, tL] <= t < 20}}]')
    eqs = M('Normal @ combined')

    return eqs

def parse_piecewise_sp(expression_string):
    # Convert Mathematica syntax to Python syntax
    expression_string = expression_string.replace('Log', 'sp.log').replace('Sqrt', 'sp.sqrt')
    expression_string = expression_string.replace("ArcTanh", "sp.atanh").replace("ArcCosh", "sp.acosh").replace("ArcSinh", "sp.asinh")
    expression_string = expression_string.replace("ArcCoth", "sp.acoth").replace("ArcSech", "sp.asech").replace("ArcCosecanth", "sp.acsch")

    # Extract individual pieces
    pieces = re.findall(r'\{(.*?),\s*Inequality\[(.*?)\]\}', expression_string, re.DOTALL)
    
    conditions = []
    expressions = []
    
    for expr, ineq in pieces:
        expr = expr.strip()
        expr = expr.replace("\n", "").replace("[", "(").replace("]", ")").replace("^", "**")
        ineq = ineq.replace("\n", "").strip()
        ineq = ineq.replace('LessEqual', '<=').replace('Less', '<').replace(',', ' and t ')  # Convert inequality format
        ineq = re.sub(r"\s*and\s*t\s*<=\s*and\s*t\s*t\s*and\s*t\s*<\s*and\s*t\s*", ",", ineq)
        ineq = ineq.split(",")
        
        conditions.append([ineq[0], ineq[1]])
        expressions.append(expr)
    expressions[0] = '0'

    if float(frac(conditions[-1][1])) < 20:
        conditions.append([conditions[-1][1], 20])
        expressions.append("0")
    
    return expressions, conditions

def graph_piece(f, cs, n_squares, index, j, i):
    #print(cs)
    # Define symbolic variable
    t = sp.Symbol('t')
    
    # Suppress warnings
    np.seterr(divide='ignore', invalid='ignore')
    
    # Generate values for plotting using `f.subs()`
    t_vals = np.linspace(0, 20, 1000)  # Define t from 0 to 20
    
    # Evaluate function at each point using `.subs()`
    y_vals = np.array([f.subs(t, val).evalf() for val in t_vals])  # Convert to numeric
    
    # Plot the function
    plt.plot(t_vals, y_vals, label=r"$PDF$", color="b")
    
    # Compute transition points
    non_diff_set = set()
    for c in cs:
        for c_ in c:
            non_diff_set.add(float(frac((c_))))  # Convert fractions if needed
    #print(non_diff_set)
    
    transition_x = np.array(list(non_diff_set))
    transition_y = np.array([f.subs(t, val-0.005).evalf() for val in transition_x])  # Compute y-values
    
    plt.scatter(transition_x, transition_y, color="red", zorder=3, label="Points", s = 10)
    plt.scatter(0, 0, color="red", zorder=3, label="Points", s = 10)
    
    # Formatting
    plt.xticks(np.arange(0, 21, 2))  # Ticks from 0 to 20 with step size of 2
    plt.xlabel("t")
    plt.title("PDF_label")
    plt.grid()
    plt.legend()
    
    # Show plot
    if j != -1:
        save_path = os.path.join("results", f"{n_squares} - {index}", f"cusp - {j}")
        os.makedirs(save_path, exist_ok=True)
        # Save the plot
        plt.savefig(os.path.join(save_path, f"graph - {i}.png"))
    else:
        save_path = os.path.join("results", f"{n_squares} - {index}", f"Final_graph.png")
        plt.savefig(save_path)
    
    # Show the plot
    plt.show()
    
    # Clear the figure to avoid overlapping plots in future calls
    plt.clf()

def integrals(x0, y0, top, bottoms, points, left):
    #M('ClearAll["Global`*"]')

    x0 = M(f"{x0}")
    y0 = M(f"{y0}")
    
    top = M(f"{top}")
    
    bottom1 = M(f"{bottoms[0]}")
    try:
        bottom2 = M(f"{bottoms[1]}")
    except:
        bottom2 = None
    try:
        bottom3 = M(f"{bottoms[2]}")
    except:
        bottom3 = None
        
    try:
        point1 = M(f"{points[0]}")
    except:
        point1 = None
    try:
        point2 = M(f"{points[1]}")
    except:
        point2 = None
        
    left = M(f"{left}")
    
    func = M('func = 1/(t*x) - x0/y0*x')
    funct = M('funct = y0/(x*(x0*x + y0*y))')
    
    if bottom2 is None:
        point1 = M('point1 = 1')
    elif bottom3 is None:
        point2 = M('point2 = 1')
    
    symbolicZeros = M('symbolicZeros = Solve[bottom1 - func == 0, x]')
    l1 = M('l1 = x /. symbolicZeros[[1]]')
    l2 = M('l2 = x /. symbolicZeros[[2]]')
    
    if bottom2 is not None:
        symbolicZeros = M('symbolicZeros = Solve[bottom2 - func == 0, x]')
        m1 = M('m1 = x /. symbolicZeros[[1]]')
        m2 = M('m2 = x /. symbolicZeros[[2]]')
    
    if bottom3 is not None:
        symbolicZeros = M('symbolicZeros = Solve[bottom3 - func == 0, x]')
        r1 = M('r1 = x /. symbolicZeros[[1]]')
        r2 = M('r2 = x /. symbolicZeros[[2]]')
    
    if bottom1 is not None and bottom2 is not None and bottom3 is not None:
        expressions = {
            "timeEnter": 'timeEnter = t /. Solve[y0 - t == 0, t][[1]]',
            "timeLeftEnd": 'timeLeftEnd = t /. Solve[l1 - left == 0, t][[1]]',
            "timeBottom1": 'timeBottom1 = t /. Solve[l2 - l1 == 0, t][[1]]',
            "timePoint1": 'timePoint1 = t /. Solve[l2 == m1, t][[1]]',
            "timeBottom2": 'timeBottom2 = t /. Solve[m2 - m1 == 0, t][[1]]',
            "timePoint2": 'timePoint2 = t /. Solve[m2 == r1, t][[1]]',
            "timeBottom3": 'timeBottom3 = t /. Solve[r2 - r1 == 0, t][[1]]',
            "timeRightEnd": 'timeRightEnd = t /. Solve[r2 - 1 == 0, t][[1]]',
            "timeLeftEndA": 'timeLeftEndA = t /. Solve[l2 - left == 0, t][[1]]',
            "timeRightEndA": 'timeRightEndA = t /. Solve[r1 - 1 == 0, t][[1]]'
        }
    
        results = {}
        for key, expr in expressions.items():
            try:
                results[key] = M(expr)
            except Exception:
                results[key] = None
            #print(f'{key}: {results[key]}')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1,
         timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA)
        
        if timeRightEnd is None and timeRightEndA is None:
            timeRightEnd = M('timeRightEnd = 20')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = 20')
        
        if timeRightEnd is None and timeRightEndA is not None:
            timeRightEnd = M('timeRightEnd = timeRightEndA')
            timeBottom3 = M('timeBottom3 = timeRightEnd')
            results['timeRightEnd'] = timeRightEndA
            results['timeBottom3'] = timeRightEndA
        if timeLeftEnd is None and timeLeftEndA is not None:
            timeLeftEnd = M('timeLeftEnd = timeLeftEndA')
            timeBottom1 = M('timeBottom1 = timeLeftEnd')
            results['timeLeftEnd'] = timeLeftEndA
            results['timeBottom1'] = timeLeftEndA

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA)
                
    elif bottom1 is not None and bottom2 is not None and bottom3 is None:
        expressions = {
            "timeEnter": 'timeEnter = t /. Solve[y0 - t == 0, t][[1]]',
            "timeLeftEnd": 'timeLeftEnd = t /. Solve[l1 - left == 0, t][[1]]',
            "timeBottom1": 'timeBottom1 = t /. Solve[l2 - l1 == 0, t][[1]]',
            "timePoint1": 'timePoint1 = t /. Solve[l2 == m1, t][[1]]',
            "timeBottom2": 'timeBottom2 = t /. Solve[m2 - m1 == 0, t][[1]]',
            "timePoint2": 'timePoint2 = t /. Solve[m2 - 1 == 0, t][[1]]',
            "timeLeftEndA": 'timeLeftEndA = t /. Solve[l2 - left == 0, t][[1]]',
            "timePoint2A": 'timePoint2A = t /. Solve[m1 - 1 == 0, t][[1]]'
        }
    
        results = {}
        for key, expr in expressions.items():
            try:
                results[key] = M(expr)
            except Exception:
                results[key] = None
            #print(f'{key}: {results[key]}')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1,
         timeBottom2, timePoint2, timeLeftEndA, timePoint2A) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeLeftEndA, timePoint2A)
        
        if timePoint2 is None and timePoint2A is None:
            timePoint2 = M('timePoint2 = 20')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = 20')
        if timePoint2 is None and timePoint2A is not None:
            timePoint2 = M('timePoint2 = timePoint2A')
            timeBottom2 = M('timeBottom2 = timePoint2A')
            results['timePoint2'] = timePoint2A
            results['timeBottom2'] = timePoint2A
        if timeLeftEnd is None and timeLeftEndA is not None:
            timeLeftEnd = M('timeLeftEnd = timeLeftEndA')
            timeBottom1 = M('timeBottom1 = timeLeftEnd')
            results['timeLeftEnd'] = timeLeftEndA
            results['timeBottom1'] = timeLeftEndA

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeLeftEndA, timePoint2A)
    
    elif bottom1 is not None and bottom2 is None and bottom3 is None:
        expressions = {
            "timeEnter": 'timeEnter = t /. Solve[y0 - t == 0, t][[1]]',
            "timeLeftEnd": 'timeLeftEnd = t /. Solve[l1 - left == 0, t][[1]]',
            "timeBottom1": 'timeBottom1 = t /. Solve[l2 - l1 == 0, t][[1]]',
            "timePoint1": 'timeRightEnd = t /. Solve[l2 - 1 == 0, t][[1]]',
            "timeLeftEndA": 'timeLeftEndA = t /. Solve[l2 - left == 0, t][[1]]',
            "timePoint1A": 'timePoint1A = t /. Solve[l1 - 1 == 0, t][[1]]'
        }
    
        results = {}
        for key, expr in expressions.items():
            try:
                results[key] = (M(expr))
            except Exception:
                results[key] = None
            #print(f'{key}: {results[key]}')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A)
        
        if timePoint1 is None and timePoint1A is None:
            timePoint1 = M('timePoint1 = 20')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = 20')

        if timePoint1 is None and timePoint1A is not None:
            timePoint1 = M('timePoint1 = timePoint1A')
            timeBottom1 = M('timeBottom1 = timePoint1')
            results['timePoint1'] = timePoint1A
            results['timeBottom1'] = timePoint1A
        if timeLeftEnd is None and timeLeftEndA is not None:
            timeLeftEnd = M('timeLeftEnd = timeLeftEndA')
            timeBottom1 = M('timeBottom1 = timeLeftEnd')
            results['timeLeftEnd'] = timeLeftEndA
            results['timeBottom1'] = timeLeftEndA

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A)
    else:
        raise ValueError("section defined improperly")
    
    try:
        f1 = M('f1 = Simplify[D[Integrate[1, {x, y0/t, 1}, {y, func, top}], t]]')
    except:
        f1 = -1
    try:
        fL = M('fL = Simplify[D[Integrate[1, {x, left, 1}, {y, func, top}], t]]')
    except:
        fL = -1
        #print("no fL")
    
    nums = set()
    for key in results.keys():
        num = results[key]
        if num is None:
            continue
        nums.add(num)
    nums = sorted(list(nums))
    print(nums)
    
    #print()
    
    function_strings = []
    function_strings.append(0)
    for num in nums:
        function_string = list('000$000$000')  # Convert string to list for mutability
        
        if timeLeftEnd <= num:
            function_string[0] = '1'
        if timeBottom1 <= num:
            function_string[1] = '1'
        if timePoint1 <= num:
            function_string[2:5] = list('1$1')  # Assign list to preserve underscore
        if bottom2 is not None:
            if timeBottom2 <= num:
                function_string[5] = '1'
            if timePoint2 <= num:
                function_string[6:9] = list('1$1')  # Assign list
        if bottom3 is not None:
            if timeBottom3 <= num:
                function_string[9] = '1'
            if timeRightEnd <= num:
                function_string[10] = '1'
    
        if bottom2 is None:
            function_string[4:] = '000$000'
        elif bottom3 is None:
            #print("true")
            function_string[8:] = '000'
        function_string = ''.join(function_string)  # Convert list back to string
        #print(function_string)  # Or store it somewhere
        function_strings.append("f$" + function_string)
    function_strings[1] = "f1"
    
    equations = []
    for function in function_strings[2:]:
        print(function)
        eq = M(f"{equations_dict[function]}")
        equations.append(eq)
    
    nums.append(0)
    nums = sorted(nums)
    combined_string = 'combined = Piecewise[{'
    for i in range(len(function_strings)):
        t0 = nums[i]
        try:
            t1 = nums[i+1]
        except:
            t1 = 20
        combined_string = combined_string + f"/{{{function_strings[i]}, {t0} <= t < {t1}}}"
        if i+1 < len(function_strings):
            combined_string = combined_string + ","
    combined_string = combined_string + "}]"
    combined_string = combined_string.replace("/", "")
    combined_string = re.sub(r"\s*--\s*", "/", combined_string)
    #print(combined_string)
    
    combined = M(f"{combined_string}")
    eqs = M('Normal @ combined')

    if bottom3 is not None:
        b1Co = M('b1Co = Integrate[funct, {x, left, point1}, {y, bottom1, top}]')
        b2Co = M('b2Co = Integrate[funct, {x, point1, point2}, {y, bottom2, top}]')
        b3Co = M('b3Co = Integrate[funct, {x, point2, 1}, {y, bottom3, top}]')
        totalVol = M('totalVol = totalVol + b1Co + b2Co + b3Co')
    elif bottom2 is not None:
        b1Co = M('b1Co = Integrate[funct, {x, left, point1}, {y, bottom1, top}]')
        b2Co = M('b2Co = Integrate[funct, {x, point1, point2}, {y, bottom2, top}]')
        totalVol = M('totalVol = totalVol + b1Co + b2Co')
    elif bottom1 is not None:
        b1Co = M('b1Co = Integrate[funct, {x, left, point1}, {y, bottom1, top}]')
        totalVol = M('totalVol = totalVol + b1Co')

    allFunctions = M('allFunctions = Append[allFunctions, eqs]')
    #print(totalVol)
        
    return eqs, totalVol