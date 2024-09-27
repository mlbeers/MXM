### RUNNING CODE
# To run the program in wurc terminal: sage MXM.py n vertical horizontal bin_width
# n should be the argument used in saddle connections: probably greater than 100,000 to get meaningful graph (int)
# vertical should be the vertical permutation as a string of the form "(1)(2,3,4)"
# horizontal the same as vertical
# bin_width should be the width of the bins you want to get the desired graph

### INITIAL SETUP
# download FileZilla
# host should be wurc.math.wisc.edu (make sure you are on school wifi or vpn)
# username is your netid
#password is your netid password
# port is 22
# this should allow you to connect to the files in your directory of the server
# you need to create a blank file called "slopes.txt" and a folder called "results" that is empty
# then you can transfer the "MXM.py" file, "slopes.txt" file and "results" by highlighting them and pressing enter
# to check the files in your directory while in terminal, type "ls"
# in terminal, you have to install pandas: sage -pip install pandas
# the code should be ready to run now

### TROUBLE SHOOTING
# if the graphs.html file displays no graphs and instead says "graph" make sure the images are in the sme directory as "MXM.py"


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

# given a sorted list of slopes, this computes the gaps between adjacent pairs and stores them in a list
def compute(seq):
    gaps = list()
    for i in range(0,len(seq) - 1):
        gap = seq[i+1] - seq[i]
        gaps.append(abs(gap))
    return gaps

# given a list of gaps and binwidth, this creates a dictionary where keys are bin gap lengths and values are how many gaps are in that bin
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

# given a sorted list of slopes and a binwidth, this function creates a scatter plot
def plot_distribution(slopes, binwidth):
    gaps = compute(slopes)
    gaps_dict = graphing_dict(gaps, binwidth)
    # this for-loop finds the smallest gap. By multiplying by 4, we find the second non-differentiable point and by multiplying by 8, we         get a good window of the distribution
    for bin in gaps_dict.keys():
        if gaps_dict[bin] != 0:
            bound = int(bin*8/binwidth)
            break
    plt.scatter(list(gaps_dict.keys())[:bound], list(gaps_dict.values())[:bound], s = 5)

# taking in arguments from terminal
num = int(sys.argv[1])
vertical = str(sys.argv[2])
horizontal = str(sys.argv[3])
binwidth = float(sys.argv[4])

# creates Sym group that is >= number of squares
Sym_group = max([len(vertical), len(horizontal)])

# converting vertical and horizontal permutation strings into raw strings to be used in re.findall
v= vertical.replace("(", "\\(")
v= v.replace(")", "\\)")
h= horizontal.replace("(", "\\(")
h= h.replace(")", "\\)")

# get all the text from the slopes file
text = open(os.path.join("results", "slopes.txt"), 'r').read()

# get the headers form the slopes.txt file
line = "v: " + str(v)+ "\nh: "+ str(h) + "\nsaddle length: " + str(num) 
replace = False

# if line is in text, that means that we have already done the computation for the slopes in this particular square-tiled surface and saddle connection number. Instead of running this long portion of code, we can scrape the text file and get our slopes list, otherwise if there is no math, we run the code to get the slopes list
if re.findall(line, text) == []:
    
    G = SymmetricGroup(Sym_group)
    s = translation_surfaces.origami(G(vertical), G(horizontal))
    
    #all saddle connections start in the bottom right of our first square, if we have a square-tiled surface that doe snot have all corners      meet at one point in a torus, then this needs to be changed (Need to find a function for this to check for such surfaces)
    sc_list = s.saddle_connections(num, 1, 0)
    
    #code directly from flatsurf that finds the slopes and gets rid of duplicates
    sc_set = set()
    for sc in sc_list:
        if sc.invert() not in sc_set:
            sc_set.add(sc)
    sc_list2 = [sc for sc in sc_set]
    slopes_all = []
    for item in sc_list2:
        slopes_all.append(item.holonomy().n())
        
    #take this list tuples and put them into a sorted list of slopes, ignoring slopes of slope zero
    slopes = set()
    for slope in slopes_all:
        if slope[0] == 0:
            continue
        x = slope[0]
        y = slope[1]
        if (x <= sqrt(num/2)) and (y <= sqrt(num/2)):
            slopes.add(y/x)
    slopes = list(slopes)
    slopes.sort()

    #write list of slopes into "slopes.txt" file with a header describing n, binwdith, vertical and horizontal permutations
    with open(os.path.join("results", "slopes.txt"), "a") as f:
       f.write("v: " + str(vertical) + "\nh: " + str(horizontal) + "\nsaddle length: " + str(num)+ "\n\n" + str(slopes) + "\n\n")

# if the slopes are already in "slopes.txt":
else:
    replace = True
    #scrape file for slopes list
    slopes = re.findall(line + r"\s+\[((?:\-?\d+\.\d+,\s)*\d+\.\d+)", text)[0]
    
    #convert this large string to a list and make floats
    slopes = slopes.split(",")
    slopes = [float(i) for i in slopes]
    slopes

#create figure
fig, ax = plt.subplots(figsize=(6, 4))
plot_distribution(slopes, binwidth)
plot_name = str(vertical) + str(horizontal) + str(num) +".png"
plt.savefig(os.path.join("results", plot_name))

#if replace is True, then the html file will replace the old graph with the new one since we changed the image tot he new graph but kept the same name
#otherwise, we need to create the html code for the header and graphs
if replace == False:
    plot_setup = """
    <html>
      <body>
        <h3> {} </h3>
        <h3> [] </h3>
        <h3> () </h3>
        <img src = plot_name alt = "graph" />
        <p>-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------</p>
      </body>
    </html>
    """
    #inserting vertical, horizonatl, saddle length, and plot name for image into the file
    plot_setup = plot_setup.replace("{}", "v: " + str(vertical))
    plot_setup = plot_setup.replace("[]", "h: " + str(horizontal))
    plot_setup = plot_setup.replace("()", "saddle length: " + str(num))
    plot_setup = plot_setup.replace("plot_name", plot_name)
    
    #write html file
    with open(os.path.join("results", "graphs.html"), "a") as f:
        f.write(plot_setup)