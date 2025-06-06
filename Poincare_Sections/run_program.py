#!/usr/bin/env python3

import argparse
import subprocess
import sys
from Library import *
import re
import os
from fractions import Fraction as frac
import dill

def run_sage(script, args, log_file):
    """Run a Sage script with the given arguments"""
    cmd = ["sage", script] + args
    print(f">>> Running: {' '.join(cmd)}")
    try:
        if log_file:
            with open(log_file, "a") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running {script}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run Sage scripts with given arguments.")

    parser.add_argument("--num_squares", type=int, help="Number of squares (positive integer greater than 5)")
    parser.add_argument("--index", type=int, help="Index (must be valid for given num_squares)")
    parser.add_argument("--horizontal_perm", type=str, help="Permutation for the horizontal pairings")
    parser.add_argument("--vertical_perm", type=str, help="Permutation for the vertical pairings")

    parser.add_argument("--dx", type=string, default='1/2000', help="dx sampling resolution (positive fraction string)")
    parser.add_argument("--dy", type=string, default='None', help="dy sampling resolution (positive fraction string)")
    parser.add_argument("--vec_length", type=int, default=2000, help="Length of the generated vectors (positive integer)")

    parser.add_argument("--folder", type=str, help="directory of the folder where you want the files to live")

    args = parser.parse_args()

    # Check mutually exclusive argument groups
    num_index_provided = args.num_squares is not None and args.index is not None
    perm_provided = args.horizontal_perm is not None and args.vertical_perm is not None

    if num_index_provided == perm_provided:
        print("Error: You must provide either --num_squares and --index OR --horizontal_perm and --vertical_perm, but not both.")
        sys.exit(1)

    # Validate num_squares and index if provided
    if num_index_provided:
        if args.num_squares <= 5:
            print("Error: --num_squares must be a greater than 5.")
            sys.exit(1)

        upper_index = len(perm_list(args.num_squares)) - 1
        if args.index < 0 or args.index > upper_index:
            print(f"Error: --index must be an index between 0 and {upper_index} for num_squares={args.num_squares}")

        try:
            perm = perms_list(args.num_squares)[args.index]
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    if num_index_provided:
        # Case 1: Using num_squares/index
        args.folder = args.folder or f"{args.num_squares}_{args.index}"  # Auto-generate if None
    else:
        # Case 2: Using permutations
        if not args.folder:
            parser.error("--folder is REQUIRED when using permutations")
        try:
            perm = Origami(args.horizontal_perm, args.vertical_perm)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    # Create folder if it doesn't exist
    os.makedirs(args.folder, exist_ok=True)  # Creates folder if needed

    # handle dx
    try:
        dx_frac = frac(args.dx)
        dx = float(dx_frac)
    except Exception as e:
        print(f"Error for dx: {e}")
        sys.exit(1)
    if dx <= 0:
        print("Error: --dx must be a positive fraction.")
        sys.exit(1)

    # handle dy
    if dy == 'None':
        dy_frac = -1
        dy = -1
    else:
        try:
            dy_frac = frac(args.dy)
            dy = float(dy_frac)
        except Exception as e:
            print(f"Error for dy: {e}")
            sys.exit(1)
        if dy <= 0:
            print("Error: --dy must be a positive fraction.")
            sys.exit(1)

    # vec_length positive int
    if args.vec_length <= 0:
        print("Error: --vec_length must be a positive integer.")
        sys.exit(1)

    # Prepare args for scripts as strings
    dx = str(args.dx)
    dy = str(args.dy)
    vec_length = str(args.vec_length)
    folder = args.folder

    # create directories if they dont exist for recording info
    os.makedirs("results", exist_ok=True)
    os.makedirs("vecs", exist_ok=True) 
    os.makedirs(os.path.join("results", folder), exist_ok=True)

    #save perm with dill to use in script
    with open(os.path.join("results", folder, "perm.dill"), 'wb') as f:
        dill.dump(perm, f)

    run_sage("script_vector.py", [vec_length, folder], os.path.join(folder, "log.txt"))
    run_sage("script_winners.py", [dx, dy, folder], os.path.join(folder, "log.txt"))
    run_sage("script_integrals.py", [folder], os.path.join(folder, "log.txt"))

    print("\nAll scripts completed successfully.")

if __name__ == "__main__":
    main()