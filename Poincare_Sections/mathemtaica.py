from sage.all import *
import numpy as np
M = mathematica
A = M('A = {{-344, 225}, {-529, 346}}')  
eigenvalues = M('eigenvalues = Eigenvalues[A]')
eigenvectors = M('eigenvectors = Eigenvectors[A]')
scaledEigenvectors = M('scaledEigenvectors = Map[LCM @@ Denominator[#]*# &, eigenvectors]')
S, J = M('{S, J} = JordanDecomposition[A]')
S_np = np.array(S.sage())

print(S_np[0][0])                                     