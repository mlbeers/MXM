# MXM
This project is part of research done at the University of Wisconsin-Madison on fine-scale statistics of straight lines on flat surfaces. It uses python code to find gap distributions for different square-tiled surfaces either explicitly or through an equivalent method using poincare sections.

## Commit instructions
Before you make changes to a fresh branch:  

`git pull` 

To create and change to a new branch (replace BRANCHNAME)  

`git checkout -b BRANCHNAME`

To change to an existing branch:  

`git checkout BRANCHNAME`

To list the current branches:  

`git branch`

At the end of the day: (make sure you're on a branch! don't work on master)  

`git status` - check which files you've changed

`git add -A` - add all the changed files to the commit

`git commit -m "your message here"` - commit the changes with a message  

`git push` - push the changes to github OR

`git push --set-upstream origin BRANCHNAME` - push a new branch to github


## Prerequisites 
- VS Code: https://code.visualstudio.com/
- Jupyter extension for VS Code: https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter
- Math IT account 


## Running through the cluster
Run this in your terminal, replacing YOURUSERNAME, and then use your Math IT credentials to log in:
`ssh -t -L 58762:localhost:58762 YOURUSERNAME@magma2.math.wisc.edu sage -n jupyter --no-browser --port-retries 0 --port=58762`

This will start a Jupyter kernel, and in the output there will be a url that looks like: http://localhost:58762... 

To run a notebook (a file that ends in .ipynb), open the file in VSCode and then in the upper right hand corner go to 
`select kernel > existing jupyter` server, paste the url from above and press enter. Now you can run the cells!

## Alternative running through cluster
ssh -t -L 64321:localhost:64321 YOURUSERNAME@wurc3.math.wisc.edu
source /usr/local/mambaforge/sf-env
jupyter notebook --no-browser --port-retries 0 --port=64321

## Running Tests
The following code will give you all the info needed about a section and shape. If you wish to test something,
copy and paste the code from script.py up until the point you wish to test:

n_squares = 7
index = 1
j = 0
dx = 0.0005
permutations = perms_list(n_squares)
perm = permutations[index]

vec_file = "vecs" + str(n_squares) + "-" + str(index) + ".npy"
vecs0 = load_arrays_from_file(os.path.join("vecs", vec_file))

with open(os.path.join("results", f"{n_squares} - {index}", "setup.dill"), 'rb') as f:
    loaded_data = dill.load(f)
a,c,e,g = loaded_data

df = read_df(n_squares, index, j)

with open(os.path.join("results", f"{n_squares} - {index}", f"secs - {j}.dill"), 'rb') as f:
    secs = dill.load(f)

## Making Test code and running scripts
Use tmux to create new windows to run code in. You can start the code in these windows and disattach so the code will run even if your computer is off
- script_vector.py generates saddle connections
- script.py does calculations with those vectors on your shape

In your window run the following:

```
cd MXM/Poincare_Sections
source /usr/local/mambaforge/sf-env
```

then run either  
- `sage script_vector n_squares index increment`
- `sage script n_squares index`

where
- n_squares is the number of squares in your STS
- index is the index of the shape you are looking at in the list of permutations
- increment is how you want to increment your index as it loops through all shapes, fine to keep as 1
