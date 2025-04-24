# Gap_Distribution - negligible

# Mathematica
## 1.nb - notebook used in Mathematica to get equations for subsections of Poincare sections with one bottom edge
## 2.nb - notebook used in Mathematica to get equations for subsections of Poincare sections with more than one bottom edge
## Functions.nb - notebook used in Mathematica to add all functions together in a list

# Poincare_Sections
## results - folder containing the results of all computations for each given STS
## vecs - folder containing the list of saddle connections for a given STS
## CoVolumne.ipynb - outdated
## Cylinders.ipynb - outdated
## Download.ipynb - notebook that downloads a shape's folder in results to a zip file that can then be transfered to the BOX
## Equations.ipynb - notebook used to compare equation equalities. Used for final equations for each STS and intermediate equations for subsections of Poincare sections.
## FirstCusp.ipynb - used to analyze the first cusp of an STS, kind of outdated
## Individual.ipynb - outdated
## Integrals.ipynb - used to test the sec2 vs sec computations to test better way of getting section dimensions
## Mathematica_Test.ipynb - used to test the Mathematica code (a bit messy)
## Optimization_Test.ipynb - outdated
## Pooling_Example.ipynb - tested how pooling works for computations in script_winners.py
## script_test.ipynb - tested different functions in script_winners.py (a bit messy)
## Testing.ipynb - testing for a lot of stuff across different functions (very messy)
## Veech.ipynb - breakdown of indexes and their veech groups across different number of squares for STSs.
## winners.ipynb - outdated
## cusp.py - outdated calculations of the sections for a given STS at a given cusp
## **integration_functions.py** - python file that contains all the "patterns" used for subsections in mathematica
- base_funcs contains all the different combinations of a parabola intersecting the bottom of a subsection at only **ONE** edge. Different edges are seperated by "$" and within each tuple of three values between "$", the first entry represents intersections with the left point of this bottom segment, second entry represents intersection with the bottom segement in general, third entry represents intersection with the right point of this bottom segment. For example, "b$000$010$000" represents a parabola intersecting the second bottom segment of an STS but not at either edge. For an STS with one bottom edge, the first tuple represents that bottom segment and the other two tuples will be 0 for all values, same logic for a shape with 2 bottom segments. **Code does not support more than three bottom segments**.
- equations_dict contains all possible combinations of base_funcs, meaning a parabola could intersect the bottom1 and bottom3 segments of a section simultaneously. Not all combinations are possible, like f$011$010$000 isnt possible, would need to be f$011$110$000 since the segments share a point, but these equations satisfy all possible equations we might have when computing gap distributions.
## **mathematica.py** - python file containing functions used in computing integrals for gap distributions. Utilizes mathematica for integral calculations.
## **Library.py** - python file containing functions used for calcualting winners for a given STS. Computes details about the sections, computes winners, determines equatiuosn that make up each subsection, get an estimated pdf for each cusp, and plot results.
## original.py - outdated
## parallelization.py - outdated
## saddle.py - outdated (I think its a dupe of vector.py)
## vector.py - python file that contains function for generated saddle connections for an STS
## utils.py - python file containing utility functions, like loading and saving arrays, not sure how much it is really used
## script.py - outdated
## **program.sh** - script that can run the following three scripts to create all the output for the gap distribution of an STS.
Args:
- num_squares - number of squares for your given STS
- index - the index of your STS in the list of permutations from the perms_list function
- dx - increment value wanted for the space between smapled points in determing winners in a Poincare sections. Increase for more complicated shapes with finer subsections. recommended value is 0.0005
## **script_vectors.py** - python file used to create saddle connection files for a given STS using fucntions from vectors.py
## **script_winners.py** - python file used to find winners and create/output Poincare sections for a given STS. Uses fucntions from Library.py. Is parallelized to create new processes that run simultaneously for each cusp of your STS
## **script_integrals.py** - python file used to produce equations for gap distributions. Uses code in mathematica.py