#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import kernel_par as ker
import potentials as pot

par = ker.params()

# Physical params
par.m         = 1 
par.Lambda    = 4.
par.C         = -100
par.hbaron2m  = 1 #197.33**2 / (2*m)

# System parm
par.E           = -60
par.Estep       = 10
par.Max_dist    = 2.
par.Precision   = 0.001
par.E_precision = 10**-7
par.l           = 0
par.nodi_voluti = 0

# discretization parm
par.N_of_bins   = 10000
par.Match       = par.N_of_bins/2

WF              = np.zeros(par.N_of_bins)     # Initiate Wave function
par.setup()


print ' ---- Numerov ----'
print ' SYSTEM: '
print ' mass    = ', par.m
print ' Cut-off = ', par.Lambda
print ' LEC     = ', par.C
print ' Hbar/2m = ', par.hbaron2m
print ' '
print ' '
  
while par.done== 0 :#and E<0:
    ker.numerov(par, WF)

print ' Energy: ', par.E
print ' -----------------'



w  = np.arange(par.Distance(0),par.Distance(par.N_of_bins), par.Step)
plt.xlabel('Distance')
plt.ylabel('Wave Function')
plt.title('Two body wave function')
plt.fill_between(w, 0, WF, facecolor='yellow', alpha=0.5)
plt.plot(w,WF)
plt.show()

