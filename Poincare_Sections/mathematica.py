from matplotlib import pyplot as plt
import os
from Library import *
from time import time
import dill
import re
import subprocess
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
def run_integrals(n_squares, index, a, perm):
    boundary_points = set()
    list_functions = []
    M('ClearAll["Global`*"]')
    totalVol = M('totalVol = 0')
    
    for j in range(len(a[0])):
        print(f"section {j}")
    
        with open(os.path.join("results", f"{n_squares}_{index}", f"secs_integrals_{j}.dill"), 'rb') as f:
            secs2 = dill.load(f)  
    
        subsections = format_subsections_for_mathematica(secs2)
    
        # Create directory for storing LaTeX outputs
        results_dir = os.path.join("results", f"{n_squares}_{index}", f"cusp_{j}")
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
                x0 = subsection["x0"]
                y0 = subsection["y0"]
                top = subsection["top"]
                bottoms = subsection["bottoms"]
                points = subsection["points"]
    
                eqs, totalVol = integrals(x0, y0, top, bottoms, points)
        
                # Convert equation to LaTeX format
                expression_string = repr(eqs[1])
                # converting from Mathematica to sympy notation
                es, cs = parse_piecewise_sp(expression_string)
                print(f"cs orig: {cs}")
                for c in cs:
                    for c_ in c:
                        boundary_points.add(c_)
    
                # creating piecewise function to be saved as a dill file, the true equation
                eq_tuples = []
                for eq, cond in zip(es, cs):
                    # converting from Mathematica to sympy notation
                    cleaned_eq = (parse_expr(eq, local_dict={'t': t, 'sp': sp}), sp.And(sp.simplify(cond[0]) <= t, t < sp.simplify(cond[1])))
                    eq_tuples.append(cleaned_eq)
    
                # save equation in sympy format from mathematica format
                f = sp.Piecewise(*eq_tuples)
                with open(os.path.join(results_dir, f"equation_{i}.dill"), 'wb') as file:
                    dill.dump(f, file)
    
                # creating piecewise function for graphing, accounting for cutting off a bound of Infinity
                eq_tuples = []
                if cs[-1][1] == "Infinity":
                    cs[-1][1] = str(max(int(float(frac(cs[-1][0])) + 10), 20))
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
    # totalVol is a global variable that is added to with every call to "integrals"
    totalVol = M('totalVol = FullSimplify[totalVol]')
    # numeric version
    totalVol2 = M('totalVolN = N[totalVol]')

    # expected coVolume
    veech_index = int(re.findall(r'index (\d+)', str(perm.veech_group()))[0])
    # numeric
    target = M(f'target = N[(Pi^2/6)*{veech_index}]')
    save_path = os.path.join("results", f"{n_squares}_{index}", "coVolume.txt")
    with open(save_path, "w") as file:
        file.write(f"{totalVol}\n\n\n")
        file.write(f"rounded coVol: {totalVol2}\n")
        file.write(f"rounded target: {target}\n")
    boundary_points = [sp.Symbol("Infinity") if bp == "Infinity" else sp.simplify(bp) for bp in boundary_points]
    sorted_numbers = sorted(boundary_points, key=lambda x: float(x) if x != sp.Symbol("Infinity") else float('inf'))

    with open(os.path.join("results", f"{n_squares}_{index}", "piecewise_list.dill"), 'wb') as f:
        dill.dump(list_functions, f)

    with open(os.path.join("results", f"{n_squares}_{index}", "boundary_points.dill"), 'wb') as f:
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
        output = {}
        output["bottoms"] = []
        output["points"] = []
        output["x0"] = "x0 = " + str(subsection.vec[0][0])
        output["y0"] = "y0 = " + str(subsection.vec[1][0])
        output["top"] = "top = " + str(subsection.top[0])
        for i in range(len(subsection.bottom)):
            output["bottoms"].append("bottom" + str(i+1) + " = " + str(subsection.bottom[i]))
        for i in range(len(subsection.points_bottom)):
            output["points"].append("point" + str(i) + " = " + str(subsection.points_bottom[i]))
        outputs.append(output)
    return outputs

# this code is not used but a clearer reference for a smple of what the integrals function computes.
# this is for subsections with only 1 bottom equation so work to decide certain things about the section is not included

# x0 is the x-component of the winning saddle
# y0 is the y-component of the winning saddle
# top is the equation of the top part of the section
# bottom1 is the equation of the bottom1 part of the section (assumes only one bottom equation for simplicity)
# left is the left-most point of the section

def nb_1(x0, y0, top, bottom1, left):
    M('ClearAll["Global`*"]')
    x0 = M(f'x0 = {x0}')
    y0 = M(f'y0 = {y0}')
    top = M(f'top = {top}')
    bottom1 = M(f'bottom1 = {bottom1}')
    left = M(f'left = {left}')
    func = M('func = 1/(t*x) - x0/y0*x')
    # the time when the parabala enters the section is always equal to y0 (always the top right point)
    timeEnter = M('timeEnter = t /. Solve[y0 - t == 0, t][[1]]')
    symbolicZeros = M('symbolicZeros = Solve[bottom1 - func == 0, x]')
    #left-most x-coordinate where parabola intersects bottom after intial intersection
    r1 = M('r1 = x /. symbolicZeros[[1]]')
    #right-most x-coordinate where parabola intersects bottom after intial intersection
    r2 = M('r2 = x /. symbolicZeros[[2]]')
    # time when the parabola hits the bottom line
    timeEdge = M(' timeEdge = t /. Solve[r2 - r1 == 0, t][[1]]')
    try:
        tR = M('tR = t /. Solve[r2 - 1 == 0, t][[1]]')
    except:
        tR = None
    
    try:
        tL = M('tL = t /. Solve[r1 - left == 0, t][[1]]')
    except:
        tL = None

    # intersection 
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
        #placefolder for infinity, used infinity to actual code
        tR = M('tR = 20')
    if tL is None and tLA is None:
        #placeholder for infinity, used infinity to actual code
        tL = M('tL = 20')
    
    if tR is None and tRA is not None:
        tR = M('tR = tRA')
        timeEdge = M('timeEdge = tR')
    if tL is None and tLA is not None:
        tL = M('tL = tLA')
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
        save_path = os.path.join("results", f"{n_squares}_{index}", f"cusp_{j}")
        os.makedirs(save_path, exist_ok=True)
        # Save the plot
        plt.savefig(os.path.join(save_path, f"graph_{i}.png"))
    else:
        save_path = os.path.join("results", f"{n_squares}_{index}", f"Final_graph.png")
        plt.savefig(save_path)
    
    # Show the plot
    plt.show()
    
    # Clear the figure to avoid overlapping plots in future calls
    plt.clf()

def integrals(x0, y0, top, bottoms, points):
    # dictionary to hold values for section
    values = {}

    values["x0"] = M(f"{x0}")
    values["y0"] = M(f"{y0}")
    
    values["top"] = M(f"{top}")

    for i in range(len(bottoms)):
        values[f"bottom{i+1}"] = M(f"{bottoms[i]}")
    for i in range(len(points)):
        values[f"point{i}"] = M(f"{points[i]}") 

    # Assumes all the sections have endpoints at x = 1
    values[f"point{len(points)}"] = M(f"point{len(points)} = 1")
        
    
    func = M('func = 1/(t*x) - x0/y0*x')
    funct = M('funct = y0/(x*(x0*x + y0*y))')

    print(values)
    expressions = {}
    expressions["timeEnter"] = 'timeEnter = t /. Solve[y0 - t == 0, t][[1]]'
    
    # find the points of intersection of each bottom segment with the parabola and solve for the times where the parabola intersects different points of the line segment
    for i in range(1, len(bottoms) + 1):
        symbolicZeros = M(f'symbolicZeros = Solve[bottom{i} - func == 0, x]')
        values[f"inter{i}$1"] = M(f'inter{i}$1 = x /. symbolicZeros[[1]]')
        values[f"inter{i}$2"] = M(f'inter{i}$2 = x /. symbolicZeros[[2]]')
        if i == 1:
            expressions["timePoint0"] = 'timePoint0 = t /. Solve[inter1$1 == point0, t][[1]]'
            expressions["timePoint0A"] = 'timePoint0A = t /. Solve[inter1$2 == point0, t][[1]]'
        if i != len(bottoms):
            expressions[f"timePoint{i}"] = f'timePoint{i} = t /. Solve[inter{i}$2 == inter{i+1}$1, t][[1]]'
            expressions[f"timePoint{i}A"] = f'timePoint{i}A = t /. Solve[inter{i}$2 == inter{i+1}$2, t][[1]]'
            expressions[f"timePoint{i}B"] = f'timePoint{i}B = t /. Solve[inter{i}$1 == inter{i+1}$1, t][[1]]'
        else:
            expressions[f"timePoint{i}"] = f'timePoint{i} = t /. Solve[inter{i}$2 == 1, t][[1]]'
            expressions[f"timePoint{i}A"] = f'timePoint{i}A = t /. Solve[inter{i}$1 == 1, t][[1]]'
        expressions[f"timeBottom{i}"] = f'timeBottom{i} = t /. Solve[inter{i}$1 - inter{i}$2 == 0, t][[1]]'

    # create these time variables in mathematica to be used later
    results = {}
    for key, expr in expressions.items():
        try:
            results[key] = (M(expr))
            print(f'{key}: {results[key]}')
        except Exception:
            results[key] = None
            print(f'{key}: None')

    # look at the 'alternate' times to account for edge cases to get correct times and document in mathematica
    for i in range(0, len(bottoms) + 1):
        if i == 0:
            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is None:
                results[f"timePoint{i}"] = M(f"timePoint{i} = Infinity")

            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is not None:
                results[f"timePoint{i}"] = M(f'timePoint{i} = timePoint{i}A')
                results[f"timeBottom{i+1}"] = M(f'timeBottom{i+1} = timePoint{i}')

        elif i == len(bottoms):
            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is None:
                results[f"timePoint{i}"] = M(f"timePoint{i} = Infinity")

            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is not None:
                results[f"timePoint{i}"] = M(f'timePoint{i} = timePoint{i}A')
                results[f"timeBottom{i}"] = M(f'timeBottom{i} = timePoint{i}')
                
        else:
            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is None and results[f"timePoint{i}B"] is None:
                results[f"timePoint{i}"] = M(f"timePoint{i} = Infinity")

            if results[f"timePoint{i}"] is None and results[f"timePoint{i}A"] is not None:
                results[f"timePoint{i}"] = M(f'timePoint{i} = timePoint{i}A')
                results[f"timeBottom{i+1}"] = M(f'timeBottom{i+1} = timePoint{i}')
        
            if results[f"timePoint{i}"] is None and results[f"timePoint{i}B"] is not None:
                results[f"timePoint{i}"] = M(f'timePoint{i} = timePoint{i}B')
                results[f"timeBottom{i}"] = M(f'timeBottom{i} = timePoint{i}')

    # base function where the parabola intersects the top portion of the section
    try:
        f1 = M('f1 = Simplify[D[Integrate[1, {x, y0/t, 1}, {y, func, top}], t]]')
    except:
        f1 = -1
    # base function where the parabola intersects the top portion OUTSIDE of the section
    try:
        fL = M('fL = Simplify[D[Integrate[1, {x, point0, 1}, {y, func, top}], t]]')
    except:
        fL = -1

    # get sorted list of critical times where the pdf will change
    nums = set()
    for key in results.keys():
        num = results[key]
        if num is None:
            continue
        nums.add(num)
    nums = sorted(list(nums))

    # functions used to subtract off excess area accumulated by f1 or fL
    # the parabola can hit the left end, interior, or right end of a line segment. f110 means the parabola has hit the left end and interior of the segment but not the right edge. f010 means the parabola has intersected the interior but neither end point
    # the extra '$i' attached to each key indexes the functions by the number of bottom segments
    basic_funcs = {}
    for i in range(1, len(bottoms)+1):
        basic_funcs[f"f000${i}"] = f'f000${i} = 0'
        basic_funcs[f"f110${i}"] = f'f110${i} = Simplify[D[Integrate[1, {{x, {f"point{i-1}"}, {f"inter{i}$2"}}}, {{y, func, {f"bottom{i}"}}}], t]]'
        basic_funcs[f"f010${i}"] = f'f010${i} = Simplify[D[Integrate[1, {{x, {f"inter{i}$1"}, {f"inter{i}$2"}}}, {{y, func, {f"bottom{i}"}}}], t]]'
        basic_funcs[f"f011${i}"] = f'f011${i} = Simplify[D[Integrate[1, {{x, {f"inter{i}$1"}, {f"point{i}"}}}, {{y, func, {f"bottom{i}"}}}], t]]'
        basic_funcs[f"f111${i}"] = f'f111${i} = Simplify[D[Integrate[1, {{x, {f"point{i-1}"}, {f"point{i}"}}}, {{y, func, {f"bottom{i}"}}}], t]]'

    # list of functions that will make up the piecewise equation for the pdf
    function_strings = []
    # if the basic_funcs key has already been called in mathematica, dont repeat the call
    visited_basic_func_string = []
    # pdf starts at zero before it enters the section
    function_strings.append(0)
    infinity = M('Infinity')
    # step through all the time intervals and figure out what critical points have been intersected and what needs to be subtracted off our base integral (f1 or fL)
    for num in nums:
        # if num is infinity, the function will be 0 so we do not need to compute a function and can exit the loop
        if num == infinity:
            break
        base = []
        # for each bottom segment, build our function string to pull the equation from the 'basic_funcs' dictionary
        for i in range(1, len(bottoms)+1):
            function_string = list(f'f000${i}')
            if results[f"timePoint{i-1}"] <= num:
                function_string[1] = '1'
            if results[f"timeBottom{i}"] <= num:
                function_string[2] = '1'
            if results[f"timePoint{i}"] <= num:
                function_string[3] = '1'  # Assign list to preserve underscore
            # determine whether or not to use f1 or fL if the left most point has already been reached for the first bottom segment (left-most)
            if i == 1:
                if function_string[1] == '1':
                    current = "fL"
                else:
                    current = "f1"
            function_string = ''.join(function_string)  # Convert list back to string
            # create variables in mathematica for each of our basic functions, skip repeats
            if function_string not in visited_basic_func_string:
                eq = M(basic_funcs[function_string])
                visited_basic_func_string.append(function_string)
            # full function string for the given time interval, includes multiple basic functions
            current = current + ' - ' + function_string
        function_strings.append(f"Simplify[{current}]")

    # add t=0 to our times
    nums.append(0)
    nums = sorted(nums)
    # full string to pass into mathematica for our piecewise equation
    combined_string = 'combined = Piecewise[{'
    print()
    print(function_strings)
    print()
    for i in range(len(nums)-1):
        t0 = nums[i]
        # make string representations of fractional times consistent
        t0 = re.sub(r"\s*-+\s*", "|", str(t0)).strip()
        # same thing but if t1 is infinity, it will throw an error
        try:
            t1 = nums[i+1]
            t1 = re.sub(r"\s*-+\s*", "|", str(t1)).strip()
        except:
            t1 = "Infinity"
        combined_string = combined_string + f"/{{{function_strings[i]}, {t0} <= t < {t1}}},"
    # remove comma at the end of string
    combined_string = combined_string[:-1]
    # add closing brackets
    combined_string = combined_string + "}]"
    # fix string representations and make fractions properly
    combined_string = combined_string.replace("/", "")
    combined_string = combined_string.replace("|", "/")
    
    combined = M(f"{combined_string}")
    eqs = M('Normal @ combined')

    # CoVolume Calculations
    # totalVol is a global variable that is not reset between function calls among sections and subsections. totalVol will be the cumalitive CoVolume for the entire STS after running this function for all cusps and will be the coVolume of the shape
    for i in range(1, len(bottoms)+1):
        totalVol = M(f'totalVol = totalVol + Integrate[funct, {{x, {f"point{i-1}"}, {f"point{i}"}}}, {{y, {f"bottom{i}"}, top}}]')

    return eqs, totalVol