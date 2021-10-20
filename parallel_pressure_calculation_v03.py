# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import math

import numpy as np
import pandas as pd
import scipy
import peakutils
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import find_peaks
from scipy.signal import find_peaks_cwt
import numpy as np
import peakutils
from scipy.signal import lfilter
from scipy.linalg import solveh_banded
from brokenaxes import brokenaxes

font = {'family' : 'Arial',
        'size'   : 22}

plt.rc('font', **font)


gamma=0.03
eta=15000
h=0.0000001
epsilon0=8.85E-12
epsilonp=2.5
A=0.0001
d=0.000001
U=40
delta = 0
#space1 = 1000
#space2 = 1000
#
#h = np.linspace(1E-6,2E-6,space1)
#q = np.linspace(3.2E6,3.5E6,space2)
space1 = 100
space2 = 100

h = np.linspace(1E-6,2E-6,space1)
q = np.linspace(0,5000000,space2)

freq = np.zeros((space1,space2))

qmax = np.zeros(space1)
qc = np.zeros(space1)
dp = np.zeros(space1)
p = np.zeros(space1)
#h in j
#q in i

NUM_COLORS = space1

cmaps = plt.get_cmap('gist_rainbow')

colours = [cmaps(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]

plt.figure(figsize=(9,9/1.618))
for i in range(space1):
    p[i] = -0.5*(epsilon0*epsilonp*(epsilonp-1)*U**2)/(epsilonp*(d-delta/2)-(epsilonp-1)*h[i])**2
    dp[i] = -(epsilon0*epsilonp*(epsilonp-1)**2*U**2)/(epsilonp*(d-delta/2)-(epsilonp-1)*h[i])**3
    qmax[i] = np.sqrt(-dp[i]/(2*gamma))
    qc[i] = np.sqrt(-dp[i]/gamma)
    
    for j in range(space2):
        freq[i,j] = -(h[i]**3/(3*eta))*(gamma*q[j]**4+dp[i]*q[j]**2)

    
    plt.plot(q,freq[i,:],color=colours[i])
    plt.xlabel('wavenumber ($m^{-1}$)')
    plt.ylabel('wave frequency ($s^{-1}$)')
    plt.xlim(0,5000000)
    plt.ylim(-100,100)
    plt.show()

fig = plt.figure(figsize=(9,9/1.618))
ax = fig.add_subplot(111)

fig2 = plt.figure(figsize=(9,9/1.618))
ax2 = fig2.add_subplot(111)

fig3 = plt.figure(figsize=(9,9/1.618))
ax3 = fig3.add_subplot(111)

for ind, c in enumerate(colours):
        
    ax.plot(h[ind],p[ind],'o',label='pressure',color=c)
    ax.set_ylabel('electrostatic pressure (Pa)')
    ax.set_xlabel('polymer height (m)')
#    ax.legend(loc='best',frameon=False)

    
    ax2.plot(h[ind],dp[ind],'o',label='dpressure',color=c)
    ax2.set_ylabel('dpressure (Pa/m)')
    ax2.set_xlabel('polymer height (m)')


    ax3.plot(h[ind],qc[ind],'o',color=c)
    ax3.plot(h[ind],qmax[ind],'o',color=c)
    ax3.set_ylabel('wavenumber ($m^{-1}$)')
    ax3.set_xlabel('polymer height (m)')

ax3.plot(h,qc,label='qc',color='k')
ax3.plot(h,qmax,ls='--',label='qmax',color='k')
ax3.legend(loc='best',frameon=False)
plt.show()
#from matplotlib import cm

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(h, q, freq, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.figure(figsize=(9,9/1.618))   
plt.contourf(h, q, freq, 1000, cmap='gist_rainbow')
plt.xlabel('polymer height (m)')
plt.ylabel('wavenumber ($m^{-1}$)')
plt.xlim(1.5e-6,2e-6)
plt.ylim(3.1e6,3.6e6)
plt.colorbar(label='wave frequency ($s^{-1}$)')
#plt.legend(loc='best',frameon=False)
#plt.show()
#plt.axis(aspect='image')




