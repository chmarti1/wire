#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f = np.linspace(.51,10,101)
x = 0.2

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

tc = 1./f
tf = tc * x
ta = tc * (1.-x)

ec = np.exp(-tc)
ef = np.exp(-tf)
ea = np.exp(-ta)

kmax = (1. - ef)/(1. - ec)
kmin = (ea - ec)/(1. - ec)
kmax_lim = x * (1 - 0.5*x/f) / (1 - 0.5/f)

ax.plot(f, kmax,'k')
ax.plot(f, kmax_lim, 'k--')
ax.plot(f, kmin, 'k')
ax.set_xlabel('Rotation Speed, $f\\tau$')
ax.set_ylabel('Temperature, $k$')
ax.set_ylim((0,1))
ax.grid('on')

fig.savefig('convection.png')
