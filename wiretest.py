# Wiretest
#   Generate pseudo data to test the wiretools pacakge

import wire
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import array,struct

class ProtoSignal(object):
    pass

class GaussianSignal(ProtoSignal):
    """GaussianSignal
    
A class for calculating wire signals due to a source with a Gaussian 
shape in x and y.

gs = GaussianSignal(x,y,sigma,A)

x,y     The x,y coordinate for the center of the gaussian shape
sigma   The standard deviation (length scale) of the signal
A       The peak of the signal (its value at x and y)

x,y, and sigma must have the same units as the domain length scale.
A has units current per unit length.

To obtain the current density in current/length local to a point,
I = gs(r,d,theta)
"""
    def __init__(self, x, y, sigma, A):
        self.x = x
        self.y = y
        self.sigma = np.broadcast_to(sigma, (2,))
        self.A = A
        
    def __call__(self, r, d, theta):
        
        r0 = d/np.cos(theta)
        dr = self.sigma/10
        rr = np.arange(r0,r+0.5*dr,dr)
        
        x = rr*np.cos(theta) - d
        y = rr*np.sin(theta)
        ii = self.A * np.exp(-(x*x + y*y)/(2*self.sigma*self.sigma))
        return np.trapz(ii,dx=dr)

class CircleSignal(ProtoSignal):
    """CircleSignal
   
A class for calculating wire signals due to a source with a uniform 
circular shape in x and y.

cs = CircleSignal(x,y,r,A)

x,y     The x,y coordinate for the center of the gaussian shape
r       The radius of the signal
A       The peak of the signal (its value in the circle)

x,y, and r must have the same units as the domain length scale.
A has units current per unit length.

To obtain the current density in current/length local to a point,
I = cs(r, d, theta)
"""
    def __init__(self, x, y, r, A):
        self.x = x
        self.y = y
        self.r = r
        self.A = A
        
    def __call__(self, r, d, theta):
        # modified horizontal distance to center
        x = d + self.x
        # Angle to center of the circle
        theta_center = np.arctan(self.y/x)
        # Distance to the circle's center
        L = np.sqrt(x*x + self.y*self.y)
        # Angle range to intersect the circle
        theta_crit = np.arcsin(self.r / L)
        # Adjust the angle relative to the centerline
        theta_adj = theta - theta_center
        if -theta_crit < theta_adj < theta_crit:
            cth = np.cos(theta_adj)
            sth = np.sin(theta_adj)
            # Calculate the minimum radius
            r0 = L * (cth - np.sqrt(cth*cth - 1 + self.r*self.r/L/L))
            dr = 2*np.sqrt(self.r*self.r - L*L*sth*sth)
            return self.A*min(max(0., r-r0), dr)
        else:
            return 0.
            

class TestSection:
    """TestSection
    
A class for calculating an integrated wire signal from a number of 
ProtoSignal elements.
"""

    def __init__(self):
        self.members = []
        
    def __call__(self, r, d, theta):
        I = 0.
        for m in self:
            I += m(r,d,theta)
        return I
        
    def __iter__(self):
        return self.members.__iter__()
        
    def __getitem__(self, index):
        return self.members[index]

    def addmember(self, member):
        """Add a member signal
    ts.addmember( member )
    
The member is a ProtoSignal instance like GaussianSignal or CircleSignal
"""
        if not isinstance(member, ProtoSignal):
            raise Exception('Members must be a ProtoSignal instance')
        self.members.append(member)
        
        
    def delmember(self, index):
        """Remove a member signal
    member = ts.delmember( index )

Remove a member from the member signal list and return it
"""
        return self.members.pop(index)

    def generate(self, filename, r, d, theta, show=False):
        """Generate a wirefile dataset
    ts.generate(filename, d, theta)
    
** filename **
The filename of the WireFile to generate

** r **
Scalar or array-like wire radii

** d, theta**
Array-like values of d and theta to use when generating the file. 
"""
        
        r = np.atleast_1d(r)
        if show:
            Nr = len(r)
            fig,ax = plt.subplots(1,Nr,squeeze=False)
            II = np.empty_like(theta, dtype=float)

        wf = wire.WireFile(filename)
        with wf.open('w'):
            for dd in d:
                for rindex, rr in enumerate(r):
                    for index,th in enumerate(theta):
                        I = self(rr,dd,th)
                        if show:
                            II[index] = I
                        wf.writeline(rr, dd, th, I)
                    if show:
                        ax[rindex,0].plot(theta,II,'k')
        if show:
            plt.show()

class WireData:
    def __init__(self, filename):
        self.N = None
        self.L = None
        self.C = None
        self.ncoef = 0
        
        with open(filename,'rb') as ff:
            raw = ff.read(struct.calcsize('II'))
            self.N = np.array(struct.unpack('II',raw), dtype=int)
            raw = ff.read(struct.calcsize('dd'))
            self.L = np.array(struct.unpack('dd', raw), dtype=float)
            self.ncoef = np.prod(2*self.N+1)
            raw = ff.read(2 * self.ncoef * struct.calcsize('d'))
            raw = array.array('d', raw)
            self.C = np.array(raw[0::2]) + 1j*np.array(raw[1::2])
    
        # Initialize the nu arrays
        self.nu = np.meshgrid(
                np.arange(-self.N[0], self.N[0]+1)/self.L[0],
                np.arange(-self.N[1], self.N[1]+1)/self.L[1])
        self.nu[0] = self.nu[0].reshape((self.ncoef,))
        self.nu[1] = self.nu[1].reshape((self.ncoef,))
    
    def __call__(self, x, y):
        x,y = np.broadcast_arrays(x,y)
        zshape = x.shape
        x,y = x.reshape((x.size,1)),y.reshape((y.size,1));
        
        nux = self.nu[0].reshape((1,self.ncoef))
        nuy = self.nu[1].reshape((1,self.ncoef))
        
        z = np.dot( np.exp(2j*np.pi*(nux*x + nuy*y)), self.C).real
        
        return z.reshape(zshape)
            
    def grid(self,Nx=None, Ny=None):
        """GRID - evaluate the solution at grid points
    x,y,I = grid()
        OR
    x,y,I = grid(N)
        Or
    x,y,I = grid(Nx, Ny)
    
If no argument is given, the grid spacing is selected automatically 
based on the highest wave number in each axis.  If a single scalar 
argument is given, it is treated as the number of grid oints in each
axis.  If tuple pair of arguments are found, they are interpreted as 
the number of grid points in the x- and y-axes.
"""
        # If no arguments are given
        if Nx is None:
            if Ny is None:
                Nx,Ny = 2*self.N+1
            else:
                raise Exception('Cannot specify Ny without specifying Nx.')
        # If Nx is given
        elif Ny is None:
            Ny = Nx
        # If Nx and Ny are given, do nothing        
        x = np.linspace(0,self.L[0],Nx)
        y = np.linspace(-self.L[1]/2, self.L[1]/2, Ny)
        x,y = np.meshgrid(x,y)
        return x,y
        

if __name__ == '__main__':
    ws = wire.WireSlice(1,k=2)
    ts = TestSection()
    ts.addmember(CircleSignal(0.5,0,0.4,1.))
    ts.addmember(CircleSignal(0.5,0,0.3,-1.))
    ts.generate('test.wf', 3, np.linspace(3,2,51), np.linspace(-.5,.5,51))

    ws.read('test.wf')
    ws.solve()
