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
from sage.all import *
from surface_dynamics.all import Origami
from surface_dynamics.all import *
import numpy as np
from fractions import Fraction as frac
M = mathematica
from IPython.display import display, Math
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from mathematica import *
from integration_functions import *
from sympy import symbols, And, Piecewise, Add, N
from sympy.core.relational import Relational
t = sp.Symbol('t')

t0 = time()

n_squares = int(sys.argv[1])
# index to start at
index = int(sys.argv[2])

perm = perms_list(n_squares)[index]
dx = 0.0005

vec_file = "vecs" + str(n_squares) + "-" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))

with open(os.path.join("results", f"{n_squares} - {index}", "setup.dill"), 'rb') as f:
    loaded_data = dill.load(f)
a,c,e,g = loaded_data

list_integrals, boundary_points = run_integrals(n_squares, index, a, perm)

# Create the combined Piecewise function
combined_pw, combined_scaled_pw = create_combined_piecewise(list_integrals, boundary_points)

# Print the result
latex_expr = sp.latex(combined_pw)
latex_expr = latex_expr.replace(r"\text{for}\: t", "").replace(r"\geq", "").replace(r"\wedge", r"\leq")

interval_list = [[boundary_points[i], boundary_points[i + 1]] for i in range(len(boundary_points) - 1)]
interval_list[-1][1] = interval_list[-1][0] + 10

graph_piece(combined_pw, interval_list, n_squares, index, -1, 50)

final_dir = os.path.join("results", f"{n_squares} - {index}")
os.makedirs(final_dir, exist_ok=True)  # Ensure directory exists

with open(os.path.join(final_dir, f"final_eq.dill"), 'wb') as file:
    dill.dump(combined_pw, file)

with open(os.path.join(final_dir, f"final_scaled_eq.dill"), 'wb') as file:
    dill.dump(combined_scaled_pw, file)

latex_file_path = os.path.join(final_dir, f"final_eq.tex")
pdf_file_path = os.path.join(final_dir, f"final_eq.pdf")

    # Write LaTeX file
with open(latex_file_path, "w") as latex_file:
    
    latex_file.write("\\documentclass{article}\n")
    latex_file.write("\\usepackage{amsmath}\n")
    latex_file.write("\\usepackage[paperheight=11in,paperwidth=60in]{geometry}\n")
    latex_file.write("\\begin{document}\n\n")
    latex_file.write(f"Equation:\n\\[\n{latex_expr}\n\\]\n\n")
    latex_file.write("\\end{document}\n")

# Compile LaTeX file
try:
    with open(os.devnull, "w") as FNULL:
        subprocess.run(
            ["pdflatex", "-output-directory", final_dir, latex_file_path],
            stdout=FNULL, stderr=FNULL, check=True
        )

    print(f"PDF successfully created: {pdf_file_path}")

    # Remove auxiliary files (.log, .aux)
    for ext in [".log", ".aux"]:
        aux_file = os.path.join(final_dir, f"final_eq{ext}")
        if os.path.exists(aux_file):
            os.remove(aux_file)

except subprocess.CalledProcessError:
    print(f"Error: Failed to compile {latex_file_path} into a PDF.")

t1 = time()
print(f"integrals done: {(t1-t0)/60**2}\n")
print("-----------------------------------------------------------------------------------")