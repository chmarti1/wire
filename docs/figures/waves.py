#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

Nx = 3
Ny = 4

x,y = np.meshgrid(np.linspace(0,1,101),np.linspace(0,1,101))


fig,ax = plt.subplots(Ny, Nx,
        sharex=True, sharey=True,
        layout='compressed')

for yi in range(Ny):
    for xi in range(Nx):
        ax[Ny-yi-1,xi].pcolor(np.sin(2*np.pi*(xi*x + yi*y)),cmap='Greys')
        ax[Ny-yi-1,xi].set_xticks([])
        ax[Ny-yi-1,xi].set_yticks([])

fig.savefig('waves.png')
