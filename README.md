# MXM
This project is part of research done at the University of Wisconsin-Madison on fine-scale statistics of straight lines on flat surfaces. It uses python code to find gap distributions for different square-tiled surfaces either explicitly or through an equivalent method using poincare sections.

# Test suite
The test suite for this code uses Python's built in unittest module. Note that 
the unittest module is not compatible with jupyter notebooks, so we use it only
for library functions in .py files. See the following link for an introductory 
reference:
https://www.digitalocean.com/community/tutorials/python-unittest-unit-test-example

A test in the test suite will test for one of these two things:
- correctness: we want to ensure that the API calls remain correct through 
  refactoring work.
- speed: this code is resource intensive, and these tests validate whether new 
  implementations are indeed faster than old ones.

These tests are called "unit tests", and they make use of "mock" data, which is 
data we save from running the code on examples we know well. See the following 
resources for more information on best practices for testing: https://microsoft.github.io/code-with-engineering-playbook/automated-testing/unit-testing/mocking/  

Some of the mock data is stored directly in the test code, and the rest is stored
in the following files:
- Poincare_Sections/vecs/test_data_4.npy
- Poincare_Sections/vecs7-4.npy

## Test suite status
Currently there are tests on the following files and functions:
#### In `Poincare_Sections/computations.py`:
- generate_permutations
- generate_vectors
- poincare_details (look for the generate_poincare_section_details wrapper)
#### In `Poincare_Sections/Poincare.py`:
- compute_winning_vecs_on_edges

## Running unit tests
run all unit tests on a file
`sage --python3 -m unittest Poincare_Sections/computations.py`

run a specific test
`sage --python3 -m unittest Poincare_Sections/computations.py -k test_generating_permutations`

This will use sage's install of python to run the python module unittest framework.

