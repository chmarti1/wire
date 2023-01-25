#!/usr/bin/python3
"""Spinning Disc Langmuir Probe (SDLP) spatial measurement utilities


"""

import numpy as np

class WireSlice(object):
    """WireSlice - ion densities in a 2D slice constructed from SDLP measurements

    ws = WireSlice(Lx, Ly, ...)
    
Defines a slice with width Lx and height Ly. No length units are presumed,
but all lengths used by WireSlice instances must be the same.  For display
purposes, a unit string may be specified using the length keyword.

    ws = WireSlice(Lx, Ly, length='mm', ...)

"""
    
    def __init__(self, Lx, Ly, N=None, k=None, length=''):
        # Force L to be a two-element vector
        self.L = np.array([Lx, Ly], dtype=float)
        self.length = length

        self.N = np.ones([2], dtype=int)   # number of wavenumbers
        self.C = None   # Coefficient matrix
        self.X = None   # Serialized coefficient vector

        if N is not None:
            if k is not None:
                raise Exception('Cannot specify both N and k simultaneously')
            self.N = 1 + np.asarray(N, dtype=int).squeeze()
            self.N = broadcast_to(self.N, [2])
        elif k is not None:
            k = np.asarray(k, dtype=float).squeeze()
            N = 1 + np.floor(self.L / np.broadcast_to(k, [2]))
            self.N = np.asarray(N, dtype=int)
        else:
            raise Exception('Either N or k must be specified.')



