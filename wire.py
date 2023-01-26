#!/usr/bin/python3
"""Spinning Disc Langmuir Probe (SDLP) spatial measurement utilities


"""

import numpy as np

class WireSlice(object):
    """WireSlice - ion densities in a 2D slice constructed from SDLP measurements

    ws = WireSlice(L, R, ...)
    
** L **
The domain size can be specified as a tuple for a rectuangular domain
    ws = WireSlice(L = (width, height), ...)
or with a scalar for a square domain
    ws = WireSlice(L = side, ...)

No length units are presumed, but all lengths used by WireSlice 
instances must be the same.  For display purposes, an optional unit 
string may be specified using the length keyword.  
    ws = WireSlice( ... length='mm', ...)
    
** R **
The wire radius from the center of disc rotation is specified by a 
scalar with the same length units used to specify the domain size, L.

** Defining the Expansion **
The number of terms used in the Fourier expansion in x (horizontal) and
y (vertical) can be specified explicitly by N or implicitly by the 
highest wavenumber, k.  Either N or k is required.

** N **
The number of non-zero wavenumbers in each axis can be independently 
specified by a tuple,
    ws = WireSlice( ... N = (Nx, Ny), ...)
or with a scalar to specify them to be the same 
    ws = WireSlice( ... N = number, ...)

The total number of terms in the expansion is
    (Nx + 1) * (Ny + 1).
One is added to each dimension to include the constant term.

** k **
The maximum wavenumber is an alternative means for determining the number
of terms in the expansion.  A wavenumber is the spatial equivalent to a
frequency.  With units 1/length, it specifies the number of wave periods
in a unit length.  Therefore, the maximum wavenumber along an axis is 
k=N/L.

The maximum wavenumber can be specified for each axis in a tuple
    ws = WireSlice( ... k = (kx, ky), ...)
or with a scalar to specify them to be the same
    ws = WireSlice( ... k = number, ...)
    
The maximum wavenumber determines the tightest spatial gradient that can
be represented in the expansion.  The inverse of the maximum wavenumber
is a measure of the spatial resolution that can be represented by the
expansion.
"""
    
    def __init__(self, L, R, N=None, k=None, length=''):
        # Force L to be a two-element vector
        
        self.L = np.broadcast_to(np.asarray(L, dtype=int), (2,))
        self.length = length
        self.R = float(R)

        self.N = np.zeros([2], dtype=int)   # number of wavenumbers

        if N is not None:
            if k is not None:
                raise Exception('Cannot specify both N and k simultaneously')
            self.N[:] = N
        elif k is not None:
            self.N[:] = np.floor(k/self.L)
        else:
            raise Exception('Either N or k must be specified.')
           
        self.size = 2*np.prod(self.N+1)

        # establish m and n vectors
        self.m = np.meshgrid(range(self.N[0]+1), range(self.N[1]+1), indexing='ij', sparse = True)
           
        # initialize a solution matrix and vectors
        self.X = None
        self.A = np.zeros((self.size, self.size), dtype=float)
        self.C = np.zeros((self.size,), dtype=float)
        
    def lam(self, d, theta):
        """LAM - calculate the lambda vector for a given d and theta
"""
        # Sines and cosines will come in handy repeatedly
        sth = np.sin(theta)
        cth = np.cos(theta)
        
        # Calculate R1
        # Since these raidii calculations can be infinite, first calculate
        # the maximum inverse of radius, which will never be zero.
        # Then, we'll invert it again.
        r1 = 1./max(1./self.R, 
                cth/(self.L[0] + d), 
                2*sth/self.L[1])
                
        # Calculate R0
        # If the wire doesn't pass into the domain, r0 should be equal 
        # to r1.
        r0 = min(r1, d/cth)
        # Finally, calculate delta-r
        dr = r1 - r0
        
        # Calculate a wavenumber matrix
        k = self.m[0]*cth/self.L[0] + self.m[1]*sth/self.L[1]
        # Force k0 to have a non-zero number.  It will be ignored later
        k[0] = -1.
        # Calculate the line integral matrix
        gam = 1. / (2*np.pi*k) * np.exp(-2j*np.pi*self.m[0]*d/self.L[0])
        gam *= (np.exp(2j*np.pi*k*r1) - np.exp(2j*np.pi*k*r0))

        # Flatten gamma to map it to the appropriate indices
        gam = gam.flatten()
        
        # OK, initialize the result
        NN = 2*(self.N[0]+1)*(self.N[1]+1)
        LAM = np.empty((NN,),dtype=float)
        LAM[0] = dr     # Constant term
        LAM[1] = 1.     # Offset term
        # Assign the real and imaginary portions to the appropriate 
        # portions of the lambda vector
        LAM[2::2] = gam.real[1:]
        LAM[3::2] = gam.imag[1:]
        
        return LAM

    def include(self, d, theta, I):
        """INCLUDE - include a datum in the model
    include(d, theta, I)
    
d   -   distance between the center of rotation and the domain edge
theta - the wire angle in radians
I   -   The total current measured in this position
"""
        lam = self.lam(d,theta)
        self.A += lam.reshape((lam.size,1)) * lam
        self.C += I * lam
        
    def solve(self):
        """SOLVE - solve the system with the data already included
"""
        self.X = np.linalg.solve(self.A, self.C)
