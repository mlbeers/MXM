#!/usr/bin/env python3

import argparse
import subprocess
import sys
from Library import *
import re
import os

def run_sage(script, args, log_file=None):
    """Run a Sage script with the given arguments, optionally logging output."""
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

    parser.add_argument("--dx", type=float, default=0.0005, help="dx sampling resolution (positive float)")
    parser.add_argument("--dy", type=float, default=None, help="dy sampling resolution (positive float)")
    parser.add_argument("--dz", type=float, default=0.01, help="sampling resolution for potential winners (positive float)")
    parser.add_argument("--vec_length", type=int, default=2000, help="Length of the generated vectors (positive integer)")

    parser.add_argument("--log", type=str, default="sage_run.log", help="Log file name (optional)")
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

    if args.dx <= 0:
        print("Error: --dx must be a positive float.")
        sys.exit(1)

    if args.dz <= 0:
        print("Error: --dz must be a positive float.")
        sys.exit(1)

    # vec_length positive int
    if args.vec_length <= 0:
        print("Error: --vec_length must be a positive integer.")
        sys.exit(1)

    # dy default if not provided
    if args.dy is None:
        args.dy = -1
    elif args.dy <= 0:
        print("Error: --dy must be positive.")
        sys.exit(1)

    # Prepare args for scripts as strings
    dx = str(args.dx)
    dy = str(args.dy)
    dz = str(args.dz)
    vec_length = str(args.vec_length)
    log_file = args.log
    folder = args.folder

    run_sage("script_vector.py", [perm, vec_length, folder], log_file)
    run_sage("script_winners.py", [perm, dx, dy, dz, folder], log_file)
    run_sage("script_integrals.py", [perm, folder], log_file)

    print("\nAll scripts completed successfully.")

if __name__ == "__main__":
    main()