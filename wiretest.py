# Wiretest
#   Generate pseudo data to test the wiretools pacakge

import wire
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

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
                        II[index] = I
                        wf.writeline(rr, dd, th, I)
                    if show:
                        ax[rindex,0].plot(theta,II,'k')
        if show:
            plt.show()

if __name__ == '__main__':
    ws = wire.WireSlice(1,k=2)
    ts = TestSection()
    ts.addmember(CircleSignal(0.5,0,0.4,1.))
    ts.addmember(CircleSignal(0.5,0,0.3,-1.))
    ts.generate('test.wf', 3, np.linspace(3,2,51), np.linspace(-.5,.5,51))

    ws.read('test.wf')
    ws.solve()
