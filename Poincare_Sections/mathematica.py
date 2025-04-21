from matplotlib import pyplot as plt
#import pandas as pd
#from flatsurf import *
import os
from Library import *
#import math
from time import time
#import copy
#from scipy import integrate
#from sympy import Symbol, solve, lambdify
#import traceback
import dill
#import sys
#import unittest
#from utils import load_arrays_from_file  # testing
import re
import subprocess
#from sage.all import *
#from surface_dynamics.all import Origami
#from surface_dynamics.all import *
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

# computations originally done with notebooks in Mathematica/ directory, now automated with python
# limitation exists if list "bottoms" has 4 or more entries, need to add/modify code to handle that instance
def run_integrals(n_squares, index, a, perm):
    boundary_points = set()
    list_functions = []
    M('ClearAll["Global`*"]')
    allFunctions = M('allFunctions= {}')
    totalVol = M('totalVol = 0')
    
    for j in range(len(a[0])):
        print(f"section {j}")
    
        with open(os.path.join("results", f"{n_squares} - {index}", f"secs_integrals - {j}.dill"), 'rb') as f:
            secs2 = dill.load(f)  
    
        subsections = format_subsections_for_mathematica(secs2)
    
        # Create directory for storing LaTeX outputs
        results_dir = os.path.join("results", f"{n_squares} - {index}", f"cusp - {j}")
        os.makedirs(results_dir, exist_ok=True)
    
        # LaTeX file path
        latex_file_path = os.path.join(results_dir, f"equations.tex")
        pdf_file_path = os.path.join(results_dir, f"equations.pdf")
    
        with open(latex_file_path, "w") as latex_file:
            # Start LaTeX document
            latex_file.write("\\documentclass{article}\n")
            latex_file.write("\\usepackage{amsmath}\n")
            latex_file.write("\\begin{document}\n\n")
    
            for i, subsection in enumerate(subsections):
                print(subsection)
                if len(subsection) == 5:
                    x0, y0, top, bottom1, left = subsection
                    bottoms = [bottom1]
                    points = []
                elif len(subsection) == 7:
                    x0, y0, top, bottom1, bottom2, point1, left = subsection
                    bottoms = [bottom1, bottom2]
                    points = [point1]
                elif len(subsection) == 9:
                    x0, y0, top, bottom1, bottom2, bottom3, point1, point2, left = subsection
                    bottoms = [bottom1, bottom2, bottom3]
                    points = [point1, point2]
    
                eqs, totalVol = integrals(x0, y0, top, bottoms, points, left)
    
                # Convert equation to LaTeX format
                expression_string = repr(eqs[1])
                # converting from Mathematica to sympy notation
                es, cs = parse_piecewise_sp(expression_string)
                print(f"cs orig: {cs}")
                for c in cs:
                    for c_ in c:
                        boundary_points.add(c_)

                # creating peicewsie function to be saved as a dill file, the true equation
                eq_tuples = []
                for eq, cond in zip(es, cs):
                    # converting from Mathematica to sympy notation
                    cleaned_eq = (parse_expr(eq, local_dict={'t': t, 'sp': sp}), sp.And(sp.simplify(cond[0]) <= t, t < sp.simplify(cond[1])))
                    eq_tuples.append(cleaned_eq)
    
                f = sp.Piecewise(*eq_tuples)
                with open(os.path.join(results_dir, f"equation_{i}.dill"), 'wb') as file:
                    dill.dump(f, file)

                # creating piecewise function for graphing, accounting for cutting off a bound of Infinity
                eq_tuples = []
                if cs[-1][1] == "Infinity":
                    cs[-1][1] = max(int(float(frac(cs[-1][0])) + 10), 20)
                    print(f"cs for graph: {cs}")
                    for eq, cond in zip(es, cs):
                        eq_tuples.append((
                            parse_expr(eq, local_dict={'t': t, 'sp': sp}),
                            sp.And(sp.simplify(cond[0]) <= t, t < sp.simplify(cond[1]))
                        ))
                    f_graph = sp.Piecewise(*eq_tuples)
                else:
                    f_graph = f
                
                list_functions.append(f)
                latex_expr = sp.latex(f)
                latex_expr = latex_expr.replace(r"\text{for}\: t", "").replace(r"\geq", "").replace(r"\wedge", r"\leq")
    
                # Write equation to LaTeX file
                latex_file.write(f"Equation {i}:\n\\[\n{latex_expr}\n\\]\n\n")
    
                # Generate graph for visualization
                graph_piece(f_graph, cs, n_squares, index, j, i)
    
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
                aux_file = os.path.join(results_dir, f"equations{ext}")
                if os.path.exists(aux_file):
                    os.remove(aux_file)
    
        except subprocess.CalledProcessError:
            print(f"Error: Failed to compile {j} into a PDF.")

    # Open the file in write mode
    # symbolic version
    totalVol = M('totalVol = FullSimplify[totalVol]')
    # numeric version
    totalVol2 = M('totalVolN = N[totalVol]')

    # expected coVolume
    veech_index = int(re.findall(r'index (\d+)', str(perm.veech_group()))[0])
    # numeric
    target = M(f'target = N[(Pi^2/6)*{veech_index}]')
    save_path = os.path.join("results", f"{n_squares} - {index}", "coVolume.txt")
    with open(save_path, "w") as file:
        file.write(f"{totalVol}\n\n\n")
        file.write(f"rounded coVol: {totalVol2}\n")
        file.write(f"rounded target: {target}\n")
    boundary_points = [sp.Symbol("Infinity") if bp == "Infinity" else sp.simplify(bp) for bp in boundary_points]
    sorted_numbers = sorted(boundary_points, key=lambda x: float(x) if x != sp.Symbol("Infinity") else float('inf'))

    with open(os.path.join("results", f"{n_squares} - {index}", "piecewise_list.dill"), 'wb') as f:
        dill.dump(list_functions, f)

    with open(os.path.join("results", f"{n_squares} - {index}", "boundary_points.dill"), 'wb') as f:
        dill.dump(sorted_numbers, f)
    
    return list_functions, sorted_numbers

def create_combined_piecewise(piecewise_list, boundary_points):
    combined_piecewise = []
    combined_scaled_piecewise = []
    # variable that keeps tracks of how many times theres a non-zero function 
    # contributing to the interval starting at 1, (usually the interval [1,2]). 
    # Used to scale the final pdf to compare across other shapes
    count = 0

    for i in range(len(boundary_points) - 1):
        lower_bound = boundary_points[i]
        upper_bound = boundary_points[i + 1]
        print(lower_bound, upper_bound)

        # Create the condition for the current interval
        interval_condition = (t >= lower_bound) & (t < upper_bound)
        if upper_bound == sp.Symbol("Infinity"):
            interval_condition = (t >= lower_bound) & (t < sp.oo)
        
        # Find all equations active in this interval
        active_equations = []
        for pw in piecewise_list:
            for expr, cond in pw.args:
                # convert Mathematica representation of infinity to sympys version for later computations
                if cond.args[1] == sp.StrictGreaterThan(sp.Symbol("Infinity"), t):
                    cond = sp.And(cond.args[0], sp.StrictLessThan(t, sp.oo))
                
                # Check if the interval is fully contained within the condition
                if sp.simplify(And(interval_condition, cond)) == interval_condition:
                    if lower_bound == 1 and expr != 0:
                        count += 1
                    active_equations.append(expr)      

        # Add the active equations together
        combined_expression = Add(*active_equations)
        combined_piecewise.append((combined_expression, interval_condition))
        # ignore interval from [0,1] so there is no division by zero error
        if count == 0:
            combined_scaled_piecewise.append((combined_expression, interval_condition))
        else:
            combined_scaled_piecewise.append((combined_expression/count, interval_condition))

    # Create the final combined Piecewise function
    return Piecewise(*combined_piecewise), Piecewise(*combined_scaled_piecewise)

# sets up a list of strings to feed into mathematica for information on the subsection with the following info:
# input is a poincare section, output is a list of string representations for each subsection
def format_subsections_for_mathematica(section):
    outputs = []
    for subsection in section:
        output = []
        output.append("x0 = " + str(subsection.vec[0][0]))
        output.append("y0 = " + str(subsection.vec[1][0]))
        output.append("top = " + str(subsection.top[0]))
        for i in range(len(subsection.bottom)):
            output.append("bottom" + str(i+1) + " = " + str(subsection.bottom[i]))
        for i in range(len(subsection.points_bottom)):
            if i == 0:
                continue
            output.append("point" + str(i) + " = " + str(subsection.points_bottom[i]))
        output.append("left = " + str(subsection.points_bottom[0]))
        outputs.append(output)
    return outputs

# this code is not used but a clearer reference for a smple of what the integrals function computes.
# this is for subsections with only 1 bottom equation so work to decide certain things about the section is not included
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
    # Convert Mathematica syntax to Python syntax for sympy
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

    # if float(frac(conditions[-1][1])) < 20:
    #     conditions.append([conditions[-1][1], 20])
    #     expressions.append("0")
    
    return expressions, conditions

def graph_piece(f, cs, n_squares, index, j, i):
    #print(cs)
    # Define symbolic variable
    t = sp.Symbol('t')
    
    # Suppress warnings
    np.seterr(divide='ignore', invalid='ignore')
    
    # Generate values for plotting using `f.subs()`
    t_vals = np.linspace(0, float(frac(cs[-1][1])), 1000)  # Define t from 0 to 20

    #print(f)
    # Evaluate function at each point using `.subs()`
    y_vals = np.array([f.subs(t, val).evalf() for val in t_vals])  # Convert to numeric
    
    # Plot the function
    plt.plot(t_vals, y_vals, label=r"$PDF$", color="b")
    
    # Compute transition points
    non_diff_set = set()
    for c in cs[:-1]:
        for c_ in c:
            non_diff_set.add(float(frac((c_))))  # Convert fractions if needed
    #print(non_diff_set)
    
    transition_x = np.array(list(non_diff_set))
    transition_y = np.array([f.subs(t, val-0.005).evalf() for val in transition_x])  # Compute y-values
    
    plt.scatter(transition_x, transition_y, color="red", zorder=3, label="Points", s = 10)
    plt.scatter(0, 0, color="red", zorder=3, label="Points", s = 10)
    
    # Formatting
    plt.xticks(np.arange(0, float(frac(cs[-1][1])), 2))  # Ticks from 0 to 20 with step size of 2
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
            "timeRightEndA": 'timeRightEndA = t /. Solve[r1 - 1 == 0, t][[1]]',
            "timePoint1A": 'timePoint1 = t /. Solve[l1 == m2, t][[1]]',
            "timePoint2A": 'timePoint2 = t /. Solve[m1 == r2, t][[1]]'
        }
    
        results = {}
        for key, expr in expressions.items():
            try:
                results[key] = (M(expr))
                print(f'{key}: {results[key]}')
            except Exception:
                results[key] = None
                print(f'{key}: None')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1,
         timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA, timePoint1A, timePoint2A) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA)
        
        if timeRightEnd is None and timeRightEndA is None:
            timeRightEnd = M('timeRightEnd = Infinity')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = Infinity')
        if timePoint1 is None and timePoint1A is None:
            timePoint1 = M('timePoint1 = Infinity')
        if timePoint2 is None and timePoint2A is None:
            timePoint2 = M('timePoint2 = Infinity')
        
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
        if timePoint1 is None and timePoint1A is not None:
            timePoint1 = M('timePoint1 = timePoint1A')
            timeBottom1 = M('timeBottom1 = timePoint1')
            results['timePoint1'] = timePoint1A
            results['timeBottom1'] = timePoint1A
        if timePoint2 is None and timePoint2A is not None:
            timePoint2 = M('timePoint2 = timePoint2A')
            timeBottom2 = M('timeBottom2 = timePoint2')
            results['timePoint2'] = timePoint2A
            results['timeBottom2'] = timePoint2A

        new_vals = [timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeBottom3, timeRightEnd, timeLeftEndA, timeRightEndA, timePoint1A, timePoint2A]
        # for val in new_vals:
        #     try:
        #         print(float(val))
        #     except:
        #         print(val)
                
    elif bottom1 is not None and bottom2 is not None and bottom3 is None:
        expressions = {
            "timeEnter": 'timeEnter = t /. Solve[y0 - t == 0, t][[1]]',
            "timeLeftEnd": 'timeLeftEnd = t /. Solve[l1 - left == 0, t][[1]]',
            "timeBottom1": 'timeBottom1 = t /. Solve[l2 - l1 == 0, t][[1]]',
            "timePoint1": 'timePoint1 = t /. Solve[l2 == m1, t][[1]]',
            "timeBottom2": 'timeBottom2 = t /. Solve[m2 - m1 == 0, t][[1]]',
            "timePoint2": 'timePoint2 = t /. Solve[m2 - 1 == 0, t][[1]]',
            "timeLeftEndA": 'timeLeftEndA = t /. Solve[l2 - left == 0, t][[1]]',
            "timePoint2A": 'timePoint2A = t /. Solve[m1 - 1 == 0, t][[1]]',
            "timePoint1A": 'timePoint1 = t /. Solve[l1 == m2, t][[1]]'
        }
    
        results = {}
        for key, expr in expressions.items():
            try:
                results[key] = (M(expr))
                print(f'{key}: {results[key]}')
            except Exception:
                results[key] = None
                print(f'{key}: None')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1,
         timeBottom2, timePoint2, timeLeftEndA, timePoint2A, timePoint1A) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeLeftEndA, timePoint2A)
        
        if timePoint2 is None and timePoint2A is None:
            timePoint2 = M('timePoint2 = Infinity')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = Infinity')
        if timePoint1 is None and timePoint1A is None:
            timePoint1 = M('timePoint1 = Infinity')
            
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
        if timePoint1 is None and timePoint1A is not None:
            timePoint1 = M('timePoint1 = timePoint1A')
            timeBottom1 = M('timeBottom1 = timePoint1')
            results['timePoint1'] = timePoint1A
            results['timeBottom1'] = timePoint1A

        new_vals = [timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeBottom2, timePoint2, timeLeftEndA, timePoint2A, timePoint1A]
        # for val in new_vals:
        #     try:
        #         #print(float(val))
        #     except:
        #         #print(val)
    
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
                print(f'{key}: {results[key]}')
            except Exception:
                results[key] = None
                print(f'{key}: None')
        
        # Unpacking the results if you still want individual variables
        (timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A) = results.values()

        #print(timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A)
        
        if timePoint1 is None and timePoint1A is None:
            timePoint1 = M('timePoint1 = Infinity')
        if timeLeftEnd is None and timeLeftEndA is None:
            timeLeftEnd = M('timeLeftEnd = Infinity')

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

        new_vals = [timeEnter, timeLeftEnd, timeBottom1, timePoint1, timeLeftEndA, timePoint1A]
        # for val in new_vals:
        #     try:
        #         #print(float(val))
        #     except:
        #         #print(val)
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
        t0 = re.sub(r"\s*-+\s*", "|", str(t0)).strip()
        try:
            t1 = nums[i+1]
            t1 = re.sub(r"\s*-+\s*", "|", str(t1)).strip()
        except:
            #attempt
            if function_strings[i] == "f$111$000$000" or function_strings[i] == "f$111$111$000" or function_strings[i] == "f$111$111$111":
                break
            t1 = "Infinity"
        combined_string = combined_string + f"/{{{function_strings[i]}, {t0} <= t < {t1}}},"
    combined_string = combined_string[:-1]
    combined_string = combined_string + "}]"
    combined_string = combined_string.replace("/", "")
    combined_string = combined_string.replace("|", "/")
    print(combined_string)
    
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