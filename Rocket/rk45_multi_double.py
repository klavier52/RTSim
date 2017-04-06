# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 02:20:21 2017

@author: riehl
"""

def rk45(params,dt,tlist,funcs,dynamics):
    time = tlist[-1]
    # Extract parameters
    etol = params['etol']
    finalT = params['finalT']
    mindt = params['mindt']
    maxdt = params['maxdt']
    # This is for the last time step to end exactly
    # In the future there can be other flags which end a simulation
    # Maybe an extra in rk45(flag)
    endflag = 0
    if finalT-time < dt:
        endflag = 1
        dt = finalT-time
    
    # For each derivative, calculate the k values according to 
    # rungga-kutta fehlberg
    k1list = k1(params,dt,time,funcs,dynamics)
    
    k2list = k2(params,dt,time,funcs,dynamics,k1list)
    k3list = k3(params,dt,time,funcs,dynamics,k1list,k2list)
    k4list = k4(params,dt,time,funcs,dynamics,k1list,k2list,k3list)
    k5list = k5(params,dt,time,funcs,dynamics,k1list,k2list,k3list,k4list)
    k6list = k6(params,dt,time,funcs,dynamics,k1list,k2list,k3list,k4list,k5list)
    
    O4list = statesO4in1(funcs,k1list,k2list,k3list,k4list,k5list,k6list)
    #print 'O4 ' + str(O4list)
    O5list = statesO5in1(funcs,k1list,k2list,k3list,k4list,k5list,k6list)
    #print 'O5 ' + str(O4list)
    Rlist = Rs(O4list,O5list,funcs,dt)
    deltamin = deltas(etol,Rlist,funcs)
    if checktol(etol,Rlist,funcs) or dt*deltamin < mindt or endflag:
        if dt > maxdt:
            dt = maxdt
        for key in funcs:
            if O4list[key] == 'temp':
                O4list[key] = getattr(dynamics, key)(params,time,O4list)
                funcs[key].append(O4list[key])
            else:
                funcs[key].append(O4list[key])
        time = time + dt
        tlist.append(time)
        if not endflag:
            dt = deltamin*dt
    else:
        dt = deltamin*dt
#    print time   
#    print dt
    return (tlist, dt, funcs)
        

def k1(params,dt,time,funcs,dynamics):
    k1list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = funcs[key[:-1]][-1]
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time,estim)
            k1list[key[:-1]] = dt*result
    return k1list

def k2(params,dt,time,funcs,dynamics,k1list):    
    k2list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = funcs[key[:-1]][-1] + (1.0/4.0)*k1list[key[:-1]]
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time + (1.0/4.0)*dt, estim)
            k2list[key[:-1]] = dt*result
    return k2list
    
def k3(params,dt,time,funcs,dynamics,k1list,k2list):
    k3list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = (funcs[key[:-1]][-1] + (3.0/32.0)*k1list[key[:-1]] 
            + (9.0/32.0)*k2list[key[:-1]])
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time + (3.0/8.0)*dt, estim)
            k3list[key[:-1]] = dt*result
    return k3list    

def k4(params,dt,time,funcs,dynamics,k1list,k2list,k3list):
    k4list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = (funcs[key[:-1]][-1] + (1932.0/2197.0)*k1list[key[:-1]]
            - (7200.0/2197.0)*k2list[key[:-1]] + (7296.0/2197.0)*k3list[key[:-1]])
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time + (12.0/13.0)*dt, estim)
            k4list[key[:-1]] = dt*result
    return k4list

def k5(params,dt,time,funcs,dynamics,k1list,k2list,k3list,k4list):
    k5list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = (funcs[key[:-1]][-1] + (439.0/216.0)*k1list[key[:-1]]
            - (8.0)*k2list[key[:-1]] + (3680.0/513.0)*k3list[key[:-1]]
            - (845.0/4101.0)*k4list[key[:-1]])
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time + dt, estim)
            k5list[key[:-1]] = dt*result
    return k5list    

def k6(params,dt,time,funcs,dynamics,k1list,k2list,k3list,k4list,k5list):
    k6list = dict.fromkeys(funcs, 'temp')
    estim = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            estim[key[:-1]] = (funcs[key[:-1]][-1] - (8.0/27.0)*k1list[key[:-1]] 
            + (2.0)*k2list[key[:-1]] - (3544.0/2565.0)*k3list[key[:-1]] 
            + (1859.0/4101.0)*k4list[key[:-1]] - (11.0/40.0)*k5list[key[:-1]])
    for key in funcs:
        if key[:-1] != '':
            result = getattr(dynamics, key)(params,time + (1.0/2.0)*dt, estim)
            k6list[key[:-1]] = dt*result
    return k6list 

def statesO4in1(funcs,k1list,k2list,k3list,k4list,k5list,k6list):
    O4list = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            result = (funcs[key[:-1]][-1] + (25.0/216.0)*k1list[key[:-1]] 
            + (1408.0/2565.0)*k3list[key[:-1]] + (2197.0/4104.0)*k4list[key[:-1]] 
            - (1.0/5.0)*k5list[key[:-1]])
            O4list[key[:-1]] = result
    return O4list

def statesO5in1(funcs,k1list,k2list,k3list,k4list,k5list,k6list):
    O5list = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            result = (funcs[key[:-1]][-1] + (16.0/135.0)*k1list[key[:-1]]
            + (6656.0/12825.0)*k3list[key[:-1]]
            + (28561.0/56430.0)*k4list[key[:-1]] - (9.0/50.0)*k5list[key[:-1]]
            + (2.0/55.0)*k6list[key[:-1]])
            O5list[key[:-1]] = result
    return O5list   

def Rs(O4list,O5list,funcs,dt):
    Rlist = dict.fromkeys(funcs, 'temp')
    for key in funcs:
        if key[:-1] != '':
            result = abs(O5list[key[:-1]] - O4list[key[:-1]])/dt
            if result == 0:
                result = 1e-16
            Rlist[key[:-1]] = result
    #        print abs(O5list[ii] - O4list[ii])/dt
    return Rlist

def deltas(etol,Rlist,funcs):
    deltalist = []
    for key in funcs:
        if key[:-1] != '':
            result = (etol/(2*Rlist[key[:-1]]))**(1.0/4.0)
            deltalist.append(result)
#        print min(deltalist)
    return min(deltalist)
    
def checktol(etol,Rlist,funcs):
    # returns 1 if pass, 0 if fail
    flag = 1
    for key in funcs:
        if key[:-1] != '':
            if Rlist[key[:-1]] > etol:
                flag = 0
                #print 'flag triggered'
                return flag
    return flag
    
    #k2 = dt*ydd(time + (1.0/4.0)*dt, state + (1.0/4.0)*k1,params)
    #k3 = dt*ydd(time + (3.0/8.0)*dt,state + (3.0/32.0)*k1 + (9.0/32.0)*k2,params)
    #k4 = dt*ydd(time + (12.0/13.0)*dt, state + (1932.0/2197.0)*k1 - 
        #(7200.0/2197.0)*k2 + (7296.0/2197.0)*k3,params)
    #k5 = dt*ydd(time + dt, state + (439.0/216.0)*k1 - (8.0)*k2 + 
        #(3680.0/513.0)*k3 - (845.0/4101.0)*k4,params)
    #k6 = dt*ydd(time + (1.0/2.0)*dt, state - (8.0/27.0)*k1 + (2.0)*k2 - 
        #(3544.0/2565.0)*k3 + (1859.0/4101.0)*k4 - (11.0/40.0)*k5,params)
    
#    statein1 = (state + (25.0/216.0)*k1 + (1408.0/2565.0)*k3 
#            + (2197.0/4104.0)*k4 - (1.0/5.0)*k5)
#    statezin1 = (state + (16.0/135.0)*k1 + (6656.0/12825.0)*k3 + 
#            (28561.0/56430.0)*k4 - (9.0/50.0)*k5 + (2.0/55.0)*k6)
    
#    R = abs(statezin1-statein1)/dt
#    delta = (etol/(2*R))**(1.0/4.0)
    
#    if R <= etol:
#        time = time + dt
#        state = statein1
#        #ii = ii + 1
#        dt = delta*dt
#    else:
#        dt = delta*dt
#    
#    return (time,state)
    
    
# Make funcs a list of strings. each string is a function