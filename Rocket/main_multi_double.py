# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 23:13:33 2017

@author: riehl
"""

import dynamics
import rk45_multi_double
import matplotlib.pyplot as plt
import numpy as np

# Simulate a ball in 1D space
params = {}
params['m'] = 20
params['g'] = 9.81
params['cd'] = 0.2
params['rho'] = 1.225
params['Az'] = np.pi*0.2**2
params['tThrust'] = [0,1,2,3,4,5,6,7,8,9,10]
params['Thrust'] = [0,100,1000,2000,2500,2500,2000,1000,100,0,0]

params['etol'] = 1e-6
etol = params['etol']
params['mindt'] = 1e-6
params['maxdt'] = 0.01
dt = 0.01
params['finalT'] = 28
finalT = params['finalT']

#states = [0,0]
t = 0
ii = 0
tlist = [t]
state1list =[]
state2list = []
state3list = []
funcs = {'zdd':['na'],'zd':[0],'z':[0]}
count = 0
print abs(tlist[-1] - finalT) > etol
while (abs(tlist[-1] - finalT) > etol):
    count = count + 1
    (tlist,dt,funcs) = rk45_multi_double.rk45(params,dt,tlist,funcs,dynamics)
    print tlist[-1]
#    if count > 200000:
#        break
#tlistfix = np.unique(tlist)
#plt.plot(tlistfix,funcs['yd'])
#plt.plot(tlistfix,funcs['y'])
#plt.plot(tlistfix[1:-1],funcs['ydd'][1:-1])
#plt.plot(funcs['x'],funcs['y'])
difft = np.diff(tlist)
plt.plot(tlist,funcs['z'])

# get rid of hard codes

# checkwithout e-16 
# cart3d- get stability derivatives
# add viscous effects
# implement in dynamics
# make problem 3D
# make problem 6D 
# verify accuracy with open rocket etc, closed form solutions etc... 

# figure out what to make as a module or class etc...
# class- geometry, dynamics, aero, thrust, gravity, integrator, monte, plotter
# terrain, weather, earth(rotrate,latlon etc..)
# avl for aero params??   
# quaternions?