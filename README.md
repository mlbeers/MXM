# MXM
This project is part of research done at the University of Wisconsin-Madison on fine-scale statistics of straight lines on flat surfaces. It uses python code to find gap distributions for different square-tiled surfaces either explicitly or through an equivalent method using poincare sections.

# Running unit tests
run all unit tests on a file
`sage --python3 -m unittest Poincare_Sections/computations.py`

run a specific test
`sage --python3 -m unittest Poincare_Sections/computations.py -k test_generating_permutations`

This will use sage's install of python to run the python module unittest framework.