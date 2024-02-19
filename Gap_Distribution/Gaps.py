import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from time import time
#!/bin/env -S sage -python
from flatsurf import *
from sage.all import *
import sys
import pandas as pd
import re
import os

def load_arrays_from_file(file_path):
    # Load arrays from the NumPy file
    arrays_list = np.load(file_path, allow_pickle=True)
    
    # Ensure each element in the list is a NumPy array
    arrays_list = [np.array(array) for array in arrays_list]
    
    return arrays_list

def compute(seq):
    gaps = list()
    for i in range(0,len(seq) - 1):
        gap = seq[i+1] - seq[i]
        gaps.append(abs(gap))
    return gaps

def graphing_dict(gaps_list, binwidth):
    #Create a dictionary
    bins_dict = {}
    #Create bins of width "binwidth" to add slope differences
    bins = list(np.arange(0, 1, binwidth))
    for bin in bins:
        bins_dict[bin] = 0
    #add slope gaps to respective bins
    for gap in gaps_list:
        for bin in bins:
            if gap < bin:
                bins_dict[bin] += 1
                break
    return bins_dict

def plot_distribution(slopes, binwidth):
    fig, ax = plt.subplots(figsize=(10, 10))
    gaps = compute(slopes)
    gaps_dict = graphing_dict(gaps, binwidth)
    for bin in gaps_dict.keys():
        if gaps_dict[bin] != 0:
            bound = int(bin*8/binwidth)
            break
    ax.scatter(list(gaps_dict.keys())[:bound], list(gaps_dict.values())[:bound], s = 5)
    plt.savefig(os.path.join("gaps", f"{n_squares} - {index}"))

def slopes(vecs):
    slopes = []
    for vec in vecs:
        if vec[0] == 0:
            continue
        slope = vec[1] / vec[0]
        slopes.append(slope)
    return slopes

n_squares = int(sys.argv[1])
index = int(sys.argv[2])
path = str(sys.argv[3])
binwidth = float(sys.argv[4])

vecs = load_arrays_from_file(os.path.join("vecs", path))
slopes = slopes(vecs)
plot_distribution(slopes, binwidth)