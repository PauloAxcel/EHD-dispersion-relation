import numpy as np
import matplotlib.pyplot as plt
import sys

gamma = 0.048
mu = 1.0
A = -1.5e-21
eps_0 = 8.85e-12
eps_1 = 2.5
eps_air = 1.0
mol = 0.001
h_0 = 100e-9
d = 200.0e-9
l_0 = 50.0e-9
V = 40.0
e = 1.6e-19
NA = 6.02e23
temp = 200+273.15
kB = 1.38e-23
k = np.sqrt(2000*e**2*NA*mol/(eps_1*eps_0*kB*temp))

size = 100
dx = 1. / size
dy = 1. / size

T = 9.0
dt = .00001
n = int(T / dt)



import sympy as sym
from sympy import symbols ,diff
from sympy.utilities.lambdify import lambdify
import math
import scipy.integrate as integrate
import scipy.special as special
from scipy import interpolate


h = symbols('h')

phi = V*(1+sym.cosh(k*h))/(1+sym.cosh(k*h)+(2*eps_air/(eps_1*k*(d-h)))*sym.sinh(k*h))

    
pot_0 = 0.5*eps_air*eps_0*(eps_air/eps_1-1)*(phi/(d-h))**2

dpot_0 = diff(pot_0,h)

ddpot_0 = diff(dpot_0,h)



#x = np.arange(0, size, 1)
#y = np.arange(0, size, 1)
x, y = np.mgrid[-1:1:20j, -1:1:20j]
z = h_0 + h_0*0.1*np.sin(x + y)

xx, yy = np.mgrid[-1:1:100j, -1:1:100j]
tck = interpolate.bisplrep(x, y, z, s=0)
h_xy = interpolate.bisplev(xx[:,0], yy[0,:], tck)

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#
#levels = MaxNLocator(nbins=15).tick_values(h_xy.min(), h_xy.max())
#
#
#fig, ax1 = plt.subplots()
#
#ini = plt.contourf(x,y,h_xy,levels=levels)
#fig.colorbar(ini, ax=ax1)
#
#plt.show()




#def show_patterns(U, ax=None):
#    ax.imshow(U, cmap=plt.cm.copper,
#              interpolation = 'bilinear',
#              extend = [-1,1,-1,1])
#    ax.set_axis_off()
    
    

step_plot = n//90

for i in range(n):

    hc = h_xy[1:-1, 1:-1]

    dhdx = np.gradient(hc , dx)
    d2hdx2 = np.gradient(dhdx[0] , dx)
    
    d3hdx2dy = np.gradient(d2hdx2[1] , dy)
    d4hdx2dy2 =  np.gradient(d3hdx2dy[1] , dy)
    
    d3hdx3 = np.gradient(d2hdx2[0] , dx)
    d4hdx4 = np.gradient(d3hdx3[0] , dx)
    
    dhdy = np.gradient(hc , dy)
    d2hdx2 = np.gradient(dhdx[0] , dx)
    
    d3hdx2dy = np.gradient(d2hdx2[1] , dy)
    d4hdx2dy2 =  np.gradient(d3hdx2dy[1] , dy)
    
    d3hdx3 = np.gradient(d2hdx2[0] , dx)
    d4hdx4 = np.gradient(d3hdx3[0] , dx)
    
    dhdy = np.gradient(hc , dy)
    d2hdy2 = np.gradient(dhdy[1] , dy)
    
    d3hdy2dx = np.gradient(d2hdy2[0] , dx)
    d4hdy2dx2 =  np.gradient(d3hdy2dx[0] , dx)
    
    d3hdy3 = np.gradient(d2hdy2[1] , dy)
    d4hdy4 = np.gradient(d3hdy3[1] , dy)
    

    
    
    func = lambdify(h,pot_0,'numpy') # returns a numpy-ready function
    pot = np.array(func(hc))
    
    func2 = lambdify(h,dpot_0,'numpy') # returns a numpy-ready function
    dpot = np.array(func2(hc))
    
    func3 = lambdify(h,ddpot_0,'numpy') # returns a numpy-ready function
    ddpot = np.array(func3(hc))
    

    dh = 1e9*dt/(3*mu)*(3*hc**2*dhdx[0]*(gamma*(d3hdx3[0] + d3hdy2dx[0])-dpot * dhdx[0])+
             hc**3*(gamma*(d4hdx4[0] + d4hdy2dx2[0]) - (ddpot * dhdx[0]**2 + dpot * d2hdx2[0]))+
             3*hc**2*dhdx[0]*(gamma*(d3hdy3[0] + d3hdx2dy[0])-dpot * dhdy[0])+
             hc**3*(gamma*(d4hdy4[0] + d4hdx2dy2[0]) - (ddpot * dhdy[0]**2 + dpot * d2hdy2[0])))
            
    
    
#    h_x = [d if math.isnan(x) else x for x in h_x]
    
#    dh = np.nan_to_num(dh)
    
    h_xy[1:-1, 1:-1] = hc + dh
    

#    h_xy[0,:] = h_xy[1,:]*(1 - sum(dh[1,:]))
#    h_xy[-1,:] = h_xy[-2,:]*(1 - sum(dh[-2,:]))
#    h_xy[:,0] = h_xy[:,1]*(1 - sum(dh[:,1]))
#    h_xy[:,-1] = h_xy[:,-2]*(1 - sum(dh[:,-2]))
    h_xy[0,:] = h_xy[1,:]
    h_xy[-1,:] = h_xy[-2,:]
    h_xy[:,0] = h_xy[:,1]
    h_xy[:,-1] = h_xy[:,-2]
    
    
        
    print(sum(sum(dh)))     
#    print(step_plot)
#    print(i)
    
        
    if i == 0 or i>160:
        levels = MaxNLocator(nbins=15).tick_values(h_xy.min(), h_xy.max())
        
        fig = plt.figure(figsize=(9,9/1.618))
        ax1 = fig.add_subplot(111, projection="3d")
        plt.title(i)
        # Plot the surface.
        surf = ax1.plot_surface(xx, yy, h_xy, cmap='nipy_spectral',
                               linewidth=0, antialiased=False)
        
        # Customize the z axis.
        ax1.set_zlim(0,d )
#        ax1.zaxis.set_major_locator(LinearLocator(10))
#        ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#        
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax1.view_init(elev=45, azim=45)
        plt.savefig('test/test'+str(i)+'.png')
        plt.show()
        
    if h_xy.max() > d:
        break
#    area = np.trapz(h_x,xx)
#    print(area)

#    if i == 99:
#        plt.plot(h_x,label='last',color='black')
#        
#plt.legend(loc='best')
#plt.show()