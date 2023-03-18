#!/usr/bin/python3
""" wiretest - Generate simulated wire data to test the wsolve binary

The TestSection class constructs wire data from Signal elements in the 
test domain.  These simple shapes of ion currents spread in space that
can be superimposed to form a more complicated ion current pattern.

The current available Signal elements are:
    CircleSignal(x,y,r,a)
    GaussianSignal(x,y,sigma,a)
See their in-line documentation for more information.

They are added to a TestSection instance using the addmember() method.
For example,
    ts = TestSection()
    ts.addmember(GaussianSignal( 0.1, 0, 0.2, 2 ))
adds a Gaussian distribution with its center at (0.1, 0), with standard
deviation 0.2, and peak ampltitude 2.  After calling the addmember()
method repeatedly to add all of the desired Signal elements, the 
TestSection instance can be queried like a function,
    I = ts(r,x,y,theta)
where
    r       := wire radius
    x,y     := disc center location
    theta   := wire angle
    I       := simulated wire current
This simulates a single measurement given a single wire location.

To generate a wire data file with a single command, see the TestSection
generate() method help.

(c)2023 Christopher Martin
"""

import wire
import numpy as np
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
        self.sigma = sigma
        self.A = A
        
    def __call__(self, r, x, y, theta):
        
        dr = self.sigma/10
        rr = np.arange(0,r+0.5*dr,dr)
        
        xx = rr*np.cos(theta) + x
        yy = rr*np.sin(theta) + y
        ii = self.A * np.exp(-(xx*xx + yy*yy)/(2*self.sigma*self.sigma))
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
        
    def __call__(self, r, x, y, theta):
        # xx,yy forms a vector to the circle center
        xx = self.x - x
        yy = self.y - y
        # Angle to center of the circle
        theta_center = np.arctan2(yy,xx)
        # Distance to the circle's center
        L = np.sqrt(xx*xx + yy*yy)
        # Angle range to intersect the circle
        theta_crit = np.arcsin(self.r / L)
        # Adjust the angle relative to the centerline
        theta_adj = theta - theta_center
        if -theta_crit < theta_adj < theta_crit:
            cth = np.cos(theta_adj)
            sth = np.sin(theta_adj)
            # Calculate the minimum radius along this path
            r0 = L * (cth - np.sqrt(cth*cth - 1 + self.r*self.r/L/L))
            # Calculate the maximum distance inside the circle along this path
            dr = 2*np.sqrt(self.r*self.r - L*L*sth*sth)
            # Select one
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
        
    def __call__(self, r, x, y, theta):
        I = 0.
        for m in self:
            I += m(r,x,y,theta)
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

    def generate(self, filename, r, x, y, theta, show=False):
        """Generate a WireData dataset
    ts.generate(filename, r, x, y, theta)
    
** filename **
The filename of the WireData to generate.  Existing files will be 
overwritten.

** r **
Scalar or array-like wire radii.  If r is not a scalar, the x, y, theta
arrays are repeated for each radius specified, as if each r-value 
specifies a separate wire mounted on a disc.  Length units are arbitrary
but must be consistent.

** x, y **
Array-like values of x, y disc center to use when generating the file. 
These arrays must be the same dimensions or broadcastable to the same 
dimensions (see Numpy's broadcasting rules).  For example, a scalar for
y and 1D array for x is permitted. See **data range** below for how data
are populated.  Length units are arbitrary, but must be consistent.

** theta **
Array-like values of wire angle (in radians) to use when generating the 
file.  

** data range **
Together, the arrays of radii, disc center location, and angles 
represent a three-dimentional test matrix that will be fully populated
in the wire data file.  It is equivalent to:
    for this_r in r:
        for this_x, this_y in zip(x,y):
            for this_theta in theta:
                ... (r,x,y,theta) code goes here ...

For example:
    ts.generate('test.wd', \
            [8, 7.9], \
            np.linspace(-9,-7,101), 
            0., 
            np.linspace(-0.18, 0.18, 101))
            
This example populates a file with 20,402 data elements representing 
data collected from two wires at radii 8cm and 7.9cm on a disc moving 
between -9cm and -7cm on the x-axis, and at angles between -0.18rad to
0.18rad.  
"""
        
        r = np.atleast_1d(r)
        if show:
            Nr = len(r)
            fig,ax = plt.subplots(1,Nr,squeeze=False)
            for rindex in range(len(r)):
                ax[0,rindex].set_xlabel('Wire Angle (rad)')
                ax[0,rindex].set_ylabel('Current (arb.)')
            II = np.empty_like(theta, dtype=float)

        x,y = np.broadcast_arrays(x,y)

        wf = wire.WireData(filename)
        with wf.open('w'):
            for xx,yy in zip(x,y):
                for rindex, rr in enumerate(r):
                    for index,th in enumerate(theta):
                        I = self(rr,xx,yy,th)
                        if show:
                            II[index] = I
                        wf.writeline(rr, xx, yy, th, I)
                    if show:
                        ax[0,rindex].plot(theta,II,'k')
        if show:
            plt.show()


if __name__ == '__main__':
    ts = TestSection()
    ts.addmember(GaussianSignal(0, 0, 0.1,  1.))
    ts.addmember(GaussianSignal(0, 0, 0.05, -1.))
    ts.generate('test.wd', 4, -np.linspace(4.4,3.4,101), 0.5, np.linspace(-.3,.0,101))

