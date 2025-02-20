#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <num_squares> <index>"
    exit 1
fi

# Assign arguments to variables
num_squares=$1
index=$2

# Call Python scripts with the arguments
sage script_vector.py "$num_squares" "$index"
sage script_winners.py "$num_squares" "$index"
sage script_integrals.py "$num_squares" "$index"

