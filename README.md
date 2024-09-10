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

