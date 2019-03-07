#
#   Wire deconvolution tools
#

import numpy as np


#=======================#
# The LineSegment class #
#=======================#

class LineSegment:
    """A 2D parametric formulation of a line based on an x,y coordinate, p0, on 
the line and a length vector dx,dy = dp.  The segment is defined by

    v(s) = p0 + dp*s
    
when p0 = (x,y), dp = (dx,dy), and s is a dimensionless location parameter
between 0  and 1.  The LineSegment is initialized by

    ls0 = LineSegment(p0=(x0,y0), dp=(dx,dy))

"""
    def __init__(self, p0, dp):
        self.p0 = np.asarray(p0)
        self.dp = np.asarray(dp)
        if self.p0.size != 2 or self.dp.size != 2:
            raise Exception('LineSegment: Expected v0 or dv to be two elements wide.')

    def __call__(self, s):
        return self.p0 + s * self.dp
        
    def intersect(self, B):
        """Calculate the intersection between self and a second line.
    sA, sB = A.intersect(B)

Parameters sA and sB are position indicators such that
    A(sA) == B(sB)
    
If 0<=sA<=1 and 0<=sB<=1, then the segments intersect.  If the two line 
segments are parallel, then sA and sB are returned as None.
"""
        a = np.zeros((2,2))
        a[:,0] = self.dv[:]
        a[:,1] = -B.dv[:]
        if np.linalg.det(
        
        b = B.v0 - self.v0
        try:
            sA,sB = np.linalg.solve(a,b)
        except:
            sA,sB = None,None
        return sA,sB
    
    
#================================#
# LineSegment Creation Functions #
#================================#
def LSrtheta(p0, R, theta):
    """Define a line segment from a starting point, a length, and an angle
    LS = LSrtheta(p0, R, theta)
    
The angle, theta, is defined relative to the x-axis in radians.
    dx = R cos(theta)
    dy = R sin(theta)
"""
    return LineSegment(p0, (R*np.cos(theta), R*np.sin(theta)))
    
def LSstop(p0, p1):
    """Define a line segment from starting and stop points.
"""
    p0 = np.asarray(p0)
    p1 = np.asarray(p1)
    return LineSegment(p0, p1 - p0)

#=================#
# Grid definition #
#=================#

class Grid:
    def __init__(self, Nx, Ny, h):
        self.Nx = Nx
        self.Ny = Ny
        self.h = h
        # if Ny is not odd
        if not (self.Ny & 0x01):
            raise Exception('Grid: Ny must be odd.')
