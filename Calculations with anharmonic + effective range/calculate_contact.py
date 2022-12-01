#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 17:29:14 2022

@author: robynlearn
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl

import scipy.io as scio
import scipy.special as scsp
import scipy.integrate as scinteg
import scipy.interpolate as scinterp
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq
from scipy.fft import rfft, rfftfreq

import os
import imageio

import BuschFunc

#plot settings
hfont = {'fontname':'Helvetica'}
fontsize = 16
fontsize_ticks = 12
fontsize_leg = 12
spine_width = 1.5


def amp_ft_integ(r,k,E):
    
    
    integ = 4*np.pi/(k)*r*BuschFunc.psi_s(r, E)*np.sin(k*r)
    #integ = (1/k)*r*psi_s(r, E)*np.sin(k*r)
    
    return integ

def high_mom(k,C):
    
    y = C/k**4
    
    return y


do_adiabatic = 1
do_psi0 = 1
do_highk = 0

#energies to calculate at
E_C = np.concatenate((np.arange(-3,0.6,0.001),np.arange(0.55,1.8,0.0002)))


#lattice depths we care about
#V_L_array = np.array([1,50,100,200,300])
V_L_array = np.array([200])

E_shift = np.ones((len(E_C),len(V_L_array)))
for i in range(len(V_L_array)):
        
    for k in range(len(E_C)):
        
        print('V_L = % .0f, %.0f' % (V_L_array[i],k))
        E_shift[k,i] = E_C[k] + BuschFunc.anharm_shift_E(E_C[k], V_L_array[i])

#adiabatic contact
if do_adiabatic == 1:
    
    C_ad = np.ones((len(E_C),len(V_L_array)))
    
    for i in range(len(V_L_array)):
    
        C_ad[:,i] = -8*np.pi*np.gradient(E_shift[:,i],1/BuschFunc.a_0_func(E_C,V_L_array[i]))
        
        C_array = np.stack((E_shift[:,i],C_ad[:,i]), axis=1)
        
        np.savetxt('C_array_lower_'+str(V_L_array[i])+'ER.csv',C_array)
        
        
if (do_psi0 == 1) or (do_highk==1):
    
    k = np.linspace(20,100,100)
    
    C_FT = np.ones(len(E_C))
    psi0 = np.ones(len(E_C))
    
    r_0 = 1e-20
    
    psi_FT = np.ones((len(k),len(E_C)))
        
    for j in range(len(E_C)):

        print(j)

        [norm_int,res_int] = np.sqrt(scinteg.quad(BuschFunc.psi_s_int, 0,50, args=(E_C[j])))
        
        if do_psi0 == 1:
            psi0[j] = (r_0*BuschFunc.psi_s(r_0,E_C[j])/norm_int)**2
            
        if do_highk == 1:
            
            for l in range(len(k)):
    
                psi_FT[l,j] = scinteg.quad(amp_ft_integ, 0, 5, args=(k[l],E_C[j]))[0] 
                psi_FT[l,j] += scinteg.quad(amp_ft_integ, 5, 10, args=(k[l],E_C[j]))[0]
                psi_FT[l,j] = psi_FT[l,j]/norm_int

    
   
            C_FT[j] = curve_fit(high_mom,k,np.abs(psi_FT[:,j])**2,p0=C_ad[j])[0]
            
fig = plt.figure(dpi=350,figsize=(6, 4))

if do_adiabatic == 1:
    for i in range(len(V_L_array)):
        plt.plot(E_shift[:,i],C_ad[:,i],'-',label=r'$V_{\mathrm{L}}$ = %.0f $E_{\mathrm{R}}$' % V_L_array[i])
        
if do_psi0 == 1:
    plt.plot(E_C,psi0*16*np.pi**2,'.',linewidth=2,label='Psi0')
    
if do_highk == 1:
    plt.plot(E_C,C_FT,marker='o',linestyle='none')
    
plt.xlabel(r'Energy ($\hbar \omega$)')
plt.ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)')

plt.legend()
        
        