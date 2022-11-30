#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:34:57 2022

@author: robynlearn
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import scipy.special as scsp
import scipy.integrate as scinteg
import scipy.interpolate as scinterp
from scipy.optimize import curve_fit

import os
import imageio

import BuschFunc

dim = 1000
n = 5

V_L = 200

E = np.linspace(0,3.7,dim) #[hbar_omega]    
a_0 = BuschFunc.a_0_func(E,V_L) #[a_h0]

for i in range(len(E)):
    E[i] += BuschFunc.anharm_shift_E(E[i], V_L)


#plot settings
hfont = {'fontname':'Helvetica'}
fontsize = 18
fontsize_ticks = 16
fontsize_leg = 12

#plot against scattering length

fig = plt.figure(dpi=350,figsize=(6, 6)) 

plt.plot(1/a_0,E,'-')
plt.plot(1/a_0, BuschFunc.E_n_array(n,dim))

plt.xlim([-30,30])
plt.ylim([0,3.5])


#plt.xticks(np.arange(-8,8.1,2),fontsize=fontsize_ticks,**hfont) 
#plt.yticks(np.arange(0,10.1,1),fontsize=fontsize_ticks,**hfont) 

plt.xlabel(r'Inverse scattering length ($a_{h}^{-1}$)',fontsize=fontsize_ticks,**hfont) 
plt.ylabel(r'Energy ($\hbar \omega$)',fontsize=fontsize_ticks,**hfont) 


dim = 500
n = 5

E = np.linspace(0,3.7,dim) #[hbar_omega]

B = BuschFunc.B_func_97(E,200)

for i in range(len(E)):
     E[i] += BuschFunc.anharm_shift_E(E[i], V_L)

fig_B = plt.figure(dpi=350,figsize=(6,4)) 

#plt.plot(B,E,linestyle='none',marker='o')
plt.plot(B,E,'-')
plt.plot(B, BuschFunc.E_n_array(n,dim))

plt.xlim([190,214])
#plt.xlim([201,205])
#plt.xlim([194,195])#plt.ylim([-2,10])
plt.ylim([0,3.7])

#plt.xticks(np.arange(190,215,4),fontsize=fontsize_ticks,**hfont) 
#plt.xticks(np.arange(202,205,1),fontsize=fontsize_ticks,**hfont) 
#plt.yticks(np.arange(-2,10.1,1),fontsize=fontsize_ticks,**hfont) 

plt.xlabel(r'B-field (G)',fontsize=fontsize_ticks,**hfont) 
plt.ylabel(r'Energy ($\hbar \omega$)',fontsize=fontsize_ticks,**hfont) 