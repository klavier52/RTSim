# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 23:23:11 2017

@author: riehl
"""
import numpy as np

def zdd(p,t,f):
    Fz = Fgz(p,t,f) + Faz(p,t,f) + Ftz(p,t,f)
    return Fz/p['m']

def zd(p,t,f):
    if f['zd'] < 0 and f['z'] < 0:
        f['zd'] = 0
    return f['zd']
    
def z(p,t,f):
    if f['z'] < 0:
        f['z'] = 0
    return f['z']
###############################################################################

def ydd(p,t,f):
    Fy = 5
    ydd = Fy/p['m']
    return ydd
    
def yd(p,t,f):
    return f['yd']
    
def y(p,t,f):
    return f['y']
    
def xdd(p,t,f):
    Fx = 5*(-f['y']) - f['xd']*abs(f['xd'])*p['cd']
    return Fx/p['m']

def xd(p,t,f):
    return f['xd']
    
def x(p,t,f):
    return f['x']
    
###############################################################################
def Fgz(p,t,f):
    Fgz = -p['m']*p['g']
    # future: make this a function of altitude, lat, lon
    return Fgz
    
def Faz(p,t,f):
    Faz = -0.5*p['rho']*f['zd']*abs(f['zd'])*p['Az']*p['cd']
    # make rho a table with altitude
    # cd based on alpha
    return Faz
    
def Ftz(p,t,f):
    Ftz = np.interp(t,p['tThrust'],p['Thrust'],0,0)#interpolate with time from thrust table
    # create a thrust table
    return Ftz
    