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
    def __init__(self, Nx, Ny, delta):
        self._Nx = int(Nx)
        self._Ny = int(Ny)
        self._delta = float(delta)
        self._yoffset = -self._delta * self._Ny / 2
    
    def node_from_xy(self, i,j):
        """Calculates the node index from the xy indices
    n = G.node_from_xy(i,j)
"""
        return self.Nx * j + i
        
    def xy_from_node(self, n):
        """Calculates the xy indices from the node index
    i,j = G.xy_from_node(n)
"""
        n = int(n)
        j = n / self.Nx
        i = n - self.Nx*j
        return i,j
        
    def node_shape(self):
        """Returns the grid dimensions 
    Nx, Ny = G.node_shape()
    
Nx and Ny are the number of nodes horizontally and vertically.
"""
        return self.Nx, self.Ny
        
    def node_size(self):
        """Returns the number of nodes
    Nn = G.node_size()

Nn is the number of nodes in the grid.
"""
        return self.Nx * self.Ny
        
    def element_size(self):
        """Returns the number of elements
    Ne = G.get_element_size()

Ne is the number of elements in the grid
"""
        return (self.Nx - 1)*(self.Ny - 1)
        
    def node_coord(self, i, j=None):
        """Return the x,y coordinates of a node
    x,y = G.node_coord(n)
        OR
    x,y = G.node_coord(i,j)
    
If only one index is given, then it is interpreted as sequential node
indexing.  If two indices are given, then they are interpreted as x-y 
indices.
"""
        if j is None:
            i,j = self.xy_from_node(i)
        return i*self._delta, self._yoffset+j*self._delta

    def a_from_rdt(self, R, d, theta):
        """Calculate the "accumulation" vector from the wire parameters
    a = G.a_from_rdt(self, R, d, theta)
    
R is the wire radius from the center of rotation
d is the distance between the center of rotation and the left-most node
theta is the wire angle from positive x in radians
"""
        # Search for the point where the wire first crosses the grid
        gridLS = LSstop(self.node_coord(0,0), self.node_coord(0,self._Ny-1))
        wireLS = LSrtheta((-d,0), R, theta)
        sw, sg = wireLS.intersect(gridLS)
        # If the segments do not intersect, just return zeros
        if sg is None or sg < 0 or sg > 1 or sw is None or sw < 0 or sw > 1:
            return np.zeros((self.node_size(),), dtype=float)
            
        # Calculate which element 
            
