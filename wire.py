#!/usr/bin/python3
"""Spinning Disc Langmuir Probe (SDLP) spatial measurement utilities

# Print help and exit
    $ wire.py -h
    # wire.py -h <command>

# Collect statistics on the wire file
    $ wire.py [-c config] stat <wirefile> <target>
    
# Construct a view of the output of wsolve
    # wire.py [-c config] view <wsolvefile> <target>
"""

import numpy as np
import array,struct
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sys
from getopt import getopt



help_text = {
    'hist':"""$ wire.py [-c config] hist <wirefile> <target>

Collect statistics on the data in the wirefile using wsolve's configuration
file.  
"""

}


class WireFile:
    def __init__(self, filename):
        self.filename = filename
        self.fd = None
        self.isread = False
        self.lineformat = '@dddd'
        self.linebytes = struct.calcsize(self.lineformat)

    def open(self, mode):
        if self.fd is not None:
            raise FileExistsError('File is already open.')
        if mode in ['r', 'rb']:
            self.fd = open(self.filename, 'rb')
            self.isread = True
        elif mode in ['w', 'wb']:
            self.fd = open(self.filename, 'wb')
            self.isread = False
        elif mode in ['a', 'ab']:
            self.fd = open(self.filename, 'ab')
            self.isread = False
        elif mode in ['x', 'xb']:
            self.fd = open(self.filename, 'xb')
            self.isread = False
        return self
    
    def close(self):
        if self.fd is None:
            raise FileExistsError('File is already closed.')
        self.fd.close()
        self.fd = None
        self.isread = False
        
    def __enter__(self):
        return self
        
    def __iter__(self):
        if not self.isread:
            raise Exception('The file is not opened in read mode.')
        return self
    
    def __next__(self):
        bb = self.fd.read(self.linebytes)
        if len(bb) < self.linebytes:
            raise StopIteration
        return struct.unpack(self.lineformat, bb)
        
    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.fd is not None:
            self.fd.close()
            self.fd = None
            self.isread = False
        return False

    def readline(self):
        """Read and return a four-element tuple
    (R, D, theta, I) = wf.readline()

R is the wire radius
D is the distance from the y-axis to the center of disc rotation
theta is the wire angle (in radians) relative to the x-axis
I is the measured wire current

Returns an empty tuple if end-of-file
"""
        if self.isread:
            bb = self.fd.read(self.linebytes)
            if len(bb) < self.linebytes:
                return ()
            return struct.unpack(self.lineformat, bb)
        raise Exception('The file is not opened in read mode.')
        
    def writeline(self, r, d, theta, I):
        """Write a single radius, diameter, angle, and current entry to the file
        
    with wf.open('w'):
        wf.writeline(r,d,theta,I)
"""
        if self.isread or self.fd is None:
            raise Exception('The file is not opened in write mode.')
        self.fd.write(struct.pack(self.lineformat, r,d,theta,I))
        
    def read(self):
        """Read all data from the file
    r,d,theta,I = wf.read()
    
Reads in numpy arrays of radius, distance, angle, and current.
"""
        if not self.isread or self.fd is None:
            raise Exception('The file is not opened in read mode.')
        R = []
        D = []
        T = []
        I = []
        for r,d,t,i in self:
            R.append(r)
            D.append(d)
            T.append(t)
            I.append(i)
        return np.array(R), np.array(D), np.array(T), np.array(I)
            
        

class SliceData:
    def __init__(self, filename):
        self.N = None
        self.L = None
        self.dshift = None
        self.C = None
        self.C_mn = None
        self.I0 = None
        self.nu = None
        self.nu_mn = None
        self.ncoef = 0
        self.filename = filename
        
        with open(filename,'rb') as ff:
            raw = ff.read(struct.calcsize('II'))
            self.N = np.array(struct.unpack('II',raw), dtype=int)
            raw = ff.read(struct.calcsize('dd'))
            self.L = np.array(struct.unpack('dd', raw), dtype=float)
            raw = ff.read(struct.calcsize('d'))
            self.dshift = struct.unpack('d', raw)[0]
            self.ncoef = np.prod(2*self.N+1)
            raw = ff.read(2 * (self.ncoef+1) * struct.calcsize('d'))
            raw = array.array('d', raw)
            nn = len(raw)
            if nn//2 != self.ncoef+1:
                raise Exception(f'SliceData: Coefficient dimension missmatch; NCOEF: {self.ncoef} NREAD: {nn//2}')
            nn -= 2
            self.C = np.array(raw[0:nn:2]) + 1j*np.array(raw[1:nn:2])
            self.C_mn = np.reshape(self.C, 2*self.N+1)
            self.I0 = raw[nn]
        
        # Initialize the nu arrays
        self.nu_mn = np.meshgrid(
                np.arange(-self.N[0], self.N[0]+1)/self.L[0],
                np.arange(-self.N[1], self.N[1]+1)/self.L[1])
        self.nu = [this.reshape((self.ncoef,)) for this in self.nu_mn]

    
    def __call__(self, x, y):
        x,y = np.broadcast_arrays(x,y)
        zshape = x.shape
        x,y = x.reshape((x.size,1)),y.reshape((y.size,1));
        
        nux = self.nu[0].reshape((1,self.ncoef))
        nuy = self.nu[1].reshape((1,self.ncoef))
        
        z = np.dot( np.exp(2j*np.pi*(nux*x + nuy*y)), self.C).real
        
        return z.reshape(zshape)
        
    def __getitem__(self, key):
        if isinstance(key,tuple):
            return self.C[self.mn_to_index(*key)]
        else:
            return self.C[key]
            
    def mn_to_index(self, m, n):
        """Calculate the serial index corresponding to m,n
    index = SD.mn_to_index(m,n)
"""        
        if abs(m) > self.N[0]:
            raise KeyError(f'SliceData[m,n] index is out of range: m={m}; Nx={self.N[0]}')
        elif abs(n) > self.N[1]:
            raise KeyError(f'SliceData[m,n] index is out of range: n={n}; Ny={self.N[1]}')
        return (m+self.N[0]) + (n+self.N[1])*(2*self.N[0]+1)

    def index_to_mn(self, index):
        """Calculate the m,n integers from the serial index
    m,n = SD.mn_to_index(index)
"""        
        if index<0 or index >= self.ncoef:
            raise KeyError(f'SliceData[index] index is out of range: index={index}; ncoef={self.ncoef}')
        n,m = divmod(index, 2*self.N[0]+1)
        return m-self.N[0], n-self.N[1]
        
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
                Nx,Ny = 4*self.N+2
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

    def show(self, x, y, ax=None):
        
        # If the axes are not specified, create a new figure and a new
        # axes.
        if ax is None:
            fig,ax = plt.subplots(1,1)
        # If the axes are specified, grab the parent figure
        else:
            fig = ax.get_figure()
        
        ax.pcolor(x,y,self(x,y), shading='auto')
        ax.set_aspect(self.L[1] / self.L[0])
        ax.set_xlim([0, self.L[0]])
        ax.set_ylim([-self.L[1]/2, self.L[1]/2])
        plt.show()

if __name__ == '__main__':
    opts,args = getopt(sys.argv[1:], 'c:h')
    configfile = 'wsolve.conf'
    config = {}
    
    # Scan the options
    for oo  in opts:
        if oo[0] == '-h':
            if len(args) == 0:
                print(__doc__)
                exit(0)
            elif args[0] in help_text:
                print(help_text[args[0]])
                exit(0) 
            else:
                print("Unrecognized command: " + args[0])
                print(__doc__)
                exit(0)
        elif oo[0] == '-c':
            configfile = oo[1]
    ##
    ## Read in the configuration file
    ##
    params = [('nthread', int), 
            ('Nx', int), ('Ny', int), 
            ('Lx', float), ('Ly', float),
            ('dshift', float)]
            
    with open(configfile,'r') as fd:
        words = fd.read().split()
    
    for pstr, ptype in params:
        pfound = words.pop(0)
        if pfound != pstr:
            raise Exception('Configuration syntax error in ' + configfile + '. Expected: ' + pstr + ' Found: ' + pfound)
        try:
            config[pstr] = ptype(words.pop(0))
            print(f'{pstr:>16s} : {config[pstr]:<}')
        except:
            raise Exception('Configuration syntax error in ' + configfile + '. Failed while parsing: ' + pstr)
            
    if args[0].lower() == 'hist':
        if len(args)!=3:
            print(help_text['hist'])
            raise Exception('After command "hist" two arguments are expected.')
            
        infile = args[1]
        target = args[2]
        # Read in the data
        with WireFile(infile).open('r') as wf:
            r,d,theta,I = wf.read()
        # Apply the configured shift
        d += config['dshift']
        # Calculate x,y coordinates of the wire tip
        x = r*np.cos(theta) - d
        y = r*np.sin(theta)
        # Calculate the histogram bin sizes
        xbin = config['Lx']/(2*config['Nx'])
        ybin = config['Ly']/(2*config['Ny'])
        # Calculate coordinate bin indices
        xi = np.array(np.floor(x / xbin), dtype=int)
        yi = np.array(np.floor(y / ybin), dtype=int)
        # Find the maximum and minimum limits
        xi_min = np.min(xi)
        xi_max = np.max(xi)
        nx = xi_max - xi_min + 1
        
        yi_min = np.min(yi)
        yi_max = np.max(yi)
        ny = yi_max - yi_min + 1
        
        count = np.zeros((ny,nx), dtype=int)
        
        for xii,yii in zip(xi,yi):
            count[yii-yi_min,xii-xi_min] += 1
            
        
        xx = np.arange(xi_min,xi_max+2)*xbin
        yy = np.arange(yi_min,yi_max+2)*ybin
        xx,yy = np.meshgrid(xx,yy)
        fig,ax = plt.subplots(1,2,figsize=(12,6),sharey=True)
        
        pc = ax[0].pcolor(xx,yy,count,shading='flat',cmap='Greys')
        ax[0].set_aspect(1.0)
        fig.colorbar(pc,ax=ax[0])
        Lx = config['Lx']
        Ly = config['Ly']/2
        ax[0].add_patch(Polygon([[0,-Ly],[0,Ly],[Lx,Ly],[Lx,-Ly]],\
                facecolor='none', edgecolor='r'))
                
        ax[1].plot(x,y,'k,')
        ax[1].add_patch(Polygon([[0,-Ly],[0,Ly],[Lx,Ly],[Lx,-Ly]],\
                facecolor='none', edgecolor='r'))
        ax[1].set_aspect(1.0)
        
        
        info = f'File: {infile}\nNx:{config["Nx"]}  Ny:{config["Ny"]}  Lx:{config["Lx"]}  Ly:{config["Ly"]}  dshift:{config["dshift"]}'
        fig.text(0.5,0.95,info,ha='center',va='center')
        
        if not target.endswith('.png'):
            target = target + '.png'
        
        fig.savefig(target)
        
    elif args[0].lower() == 'view':
        infile = args[1]
        
    else:
        raise Exception('Unrecognized command: ' + args[0])
