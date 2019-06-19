# Wiretest
#   Generate pseudo data to test the wiretools pacakge

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import scipy.sparse.linalg as linalg
import wiretools as wt


Nx = 11     # 
Ny = 11     # 
delta = .1  # grid size
R = 2.      # Wire radius
Dw = .01    # Wire diameter

# Set up a grid
G = wt.Grid(Nx,Ny,delta)

if False:
    print "Generating X"
    # Initialize an empty solution matrix
    X = sparse.lil_matrix((G.size(),1), dtype=float)
    # Create a square patch in the center of the domain
    for ii in range(3,8):
        for jj in range(3,8):
            X[G.ij_to_n(ii,jj)] = 1.
            
    X = sparse.csr_matrix(X)

if False:
    print "Generating pseudo data"
    # Start generating pseudo data
    A = sparse.lil_matrix((G.size(),G.size()), dtype=float)
    B = sparse.lil_matrix((G.size(),1), dtype=float)

    for dd in np.arange(.05, 2.05, .05):
        for theta in np.arange(-np.pi/2, np.pi/2, np.pi/100):
            L = G.lam(R, dd, theta)
            I = np.pi*Dw*(L.T*X)[0,0]
            B += I * L
            A += np.pi*Dw*(L * L.T)

    A = sparse.csr_matrix(A)
    B = sparse.csr_matrix(B)

if False:
    print "Solving"
    X2 = linalg.spsolve(A,B)

# Reformat the data for visualization
XX = np.ndarray(G.N, dtype=float)
XX2 = np.ndarray(G.N, dtype=float)

for ii in range(G.N[0]):
    for jj in range(G.N[1]):
        XX[ii,jj] = X[G.ij_to_n(ii,jj), 0]
        XX2[ii,jj] = X2[G.ij_to_n(ii,jj)]
        
f = plt.figure(1)
f.clf()
plt.pcolor(XX)

f = plt.figure(2)
f.clf()
plt.pcolor(XX2)

plt.show(block=False)
