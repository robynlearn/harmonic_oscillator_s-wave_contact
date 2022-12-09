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

#constants
hbar = 6.62607015e-34/(2*np.pi)
m = 39.96399848*1.66053906660e-27

lambda_L = 1054e-9 #laser wavelength [m]
k_L = 2*np.pi/lambda_L
E_R = hbar**2*k_L**2/(2*m)
a_L = lambda_L/2


V_L = 60
omega = 2*E_R*np.sqrt(V_L)/hbar # [Hz]
a_ho = np.sqrt(2*hbar/(m*omega)) #[m]  

E = np.linspace(-10,3.7,dim) #[hbar_omega]    
#a_0 = BuschFunc.a_0_func(E,V_L) #[a_h0]
a_0 = BuschFunc.a_E_func(E)

shift = np.ones(len(E))
E_shift = np.ones(len(E))
#"""
for i in range(len(E)):
    
    shift[i] = -1*V_L*(E_R/(hbar*omega))*(a_ho*k_L)**4*BuschFunc.anharm_shift(E[i])
   # E_shift[i] = E[i] + BuschFunc.anharm_shift_E(E[i], V_L)"""
#shift = -3/2*E_R/(hbar*omega)
E_shift = E + shift

fig = plt.figure(dpi=350,figsize=(6, 6))
plt.plot(E,-1*shift*hbar*omega/E_R)
plt.xlim([1.5,3.5])
plt.ylim([1.3,4.2])

plt.xlabel(r'Energy ($\hbar \omega$)')
plt.ylabel(r'Anharmonic shift ($-E_{\mathrm{R}}$)')


#plot settings
hfont = {'fontname':'Helvetica'}
fontsize = 18
fontsize_ticks = 16
fontsize_leg = 12

#plot against scattering length

fig = plt.figure(dpi=350,figsize=(6, 6)) 

plt.plot(1/a_0,E_shift,'-')
plt.plot(1/a_0, BuschFunc.E_n_array(n,dim))

plt.xlim([-30,30])
plt.ylim([-10,3.5])


#plt.xticks(np.arange(-8,8.1,2),fontsize=fontsize_ticks,**hfont) 
#plt.yticks(np.arange(0,10.1,1),fontsize=fontsize_ticks,**hfont) 

plt.xlabel(r'Inverse scattering length ($a_{h}^{-1}$)',fontsize=fontsize_ticks,**hfont) 
plt.ylabel(r'Energy ($\hbar \omega$)',fontsize=fontsize_ticks,**hfont) 




fig_B, ax = plt.subplots(dpi=350,figsize=(4,4)) 

B = BuschFunc.B_func_97(a_0,V_L)

B_ind = np.where(((B[0:-1] - B[1:]) > 0 ))[0]

plt.plot(B,BuschFunc.E_n_array(n,dim))

#plt.plot(B,E,'c-')
#plt.plot(B,E_shift,'c--')

for ind in range(len(B_ind)):
    
    if ind == 0:

        #plt.plot(B,E,linestyle='none',marker='o')
        plt.plot(B[:B_ind[ind]],E[:B_ind[ind]],'c-')
        plt.plot(B[:B_ind[ind]],E_shift[:B_ind[ind]],'c--')
        
    if ind+1 == len(B_ind):
        
        #plt.plot(B,E,linestyle='none',marker='o')
        plt.plot(B[B_ind[ind-1]+1:B_ind[ind]],E[B_ind[ind-1]+1:B_ind[ind]],'c-')
        plt.plot(B[B_ind[ind-1]+1:B_ind[ind]],E_shift[B_ind[ind-1]+1:B_ind[ind]],'c--')
        
        #plt.plot(B,E,linestyle='none',marker='o')
        plt.plot(B[B_ind[ind]+1:],E[B_ind[ind]+1:],'c-')
        plt.plot(B[B_ind[ind]+1:],E_shift[B_ind[ind]+1:],'c--')
        
    else:

        #plt.plot(B,E,linestyle='none',marker='o')
        plt.plot(B[B_ind[ind-1]+1:B_ind[ind]],E[B_ind[ind-1]+1:B_ind[ind]],'c-')
        plt.plot(B[B_ind[ind-1]+1:B_ind[ind]],E_shift[B_ind[ind-1]+1:B_ind[ind]],'c--')        
        
       
        
plt.xlim([180,220])
#plt.xlim([201,205])
#plt.xlim([194,195])#plt.ylim([-2,10])
plt.ylim([0,3.7])

#plt.xticks(np.arange(190,215,4),fontsize=fontsize_ticks,**hfont) 
#plt.xticks(np.arange(202,205,1),fontsize=fontsize_ticks,**hfont) 
plt.yticks(np.arange(0.1,3.8,0.2),fontsize=fontsize_ticks,**hfont) 
ax.yaxis.set_ticks_position('both')
ax.tick_params(axis='y',direction='in')
plt.grid()

plt.xlabel(r'B-field (G)',fontsize=fontsize_ticks,**hfont) 
plt.ylabel(r'Energy ($\hbar \omega$)',fontsize=fontsize_ticks,**hfont) 