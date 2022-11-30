#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 17:44:37 2022

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
        
    
      


#E_C = np.linspace(1.5+0.000001,3.5-0.0000001,50)
#E_C = np.linspace(1.5+0.01,3.5-0.01,50)

#E_C = np.linspace(-10,1.5-0.01,50)
E_C = np.concatenate((np.arange(-3,0.6,0.1),np.arange(0.55,1.66,0.05)))
k = np.linspace(20,100,100)
#k = np.linspace(100,250,100)


#adiabatic contact
dEda = np.gradient(E_C,1/BuschFunc.a_0_func(E_C))

dEdB = np.gradient(E_C,BuschFunc.B_func(BuschFunc.a_0_func(E_C),200))

C_ad = -8*np.pi*dEda



#high momentum tail + short-range wf
fig_mom = plt.figure(dpi=350, figsize=(8,4))
gs = gridspec.GridSpec(1,2)
ax0 = plt.subplot(gs[0])

E_one = 0.5

T = 0.01
N = 1000
r = np.linspace(0.00001,N*T,N,endpoint=False)

k_one = rfftfreq(N,T)[:N//2]

[norm_int,res_int] = np.sqrt(scinteg.quad(BuschFunc.psi_s_int, 0,50, args=(E_one)))
#psi_FT_one = np.ones(len(k_one))
"""
for ik in range(len(k_one)):
    psi_FT_one[ik] = scinteg.quad(amp_ft_integ, 0, 10, args=(k_one[ik],E_one))[0]/norm_int
    
  """
psi_FT_one = 2*np.pi*rfft(r**2*BuschFunc.psi_s(r,E_one)/norm_int)[:N//2] 

#popt, pcov = curve_fit(high_mom,k_one[200:],(psi_FT_one[200:])**2)
popt, pcov = curve_fit(high_mom,k_one[100:400],(psi_FT_one[100:400])**2)


ax0.plot(k_one,np.abs(psi_FT_one[:N//2])**2)
ax0.plot(k_one,high_mom(k_one,*popt),'--')

k_one = np.linspace(0,100,1000)
psi_FT_one = np.ones(len(k_one))

for ik in range(len(k_one)):
    psi_FT_one[ik] = scinteg.quad(amp_ft_integ, 0, 10, args=(k_one[ik],E_one))[0]/norm_int
    
popt, pcov = curve_fit(high_mom,k_one[200:],(psi_FT_one[200:])**2)

ax0.plot(k_one,np.abs(psi_FT_one)**2)
ax0.plot(k_one,high_mom(k_one,*popt),'--')


plt.xscale('log')
plt.yscale('log')
#plt.ylim([1e-7,100])
plt.xlim([1e-1,100])
plt.xticks(fontsize=fontsize_ticks,**hfont)
plt.yticks(fontsize=fontsize_ticks,**hfont)

plt.xlabel(r'$k$ $(a_{\mathrm{ho}}^{-1})$',fontsize=fontsize,**hfont)
plt.ylabel(r'$n(k)$ (arb.)',fontsize=fontsize,**hfont)

plt.setp(ax0.spines.values(), linewidth=spine_width)
ax0.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)

#plt.plot(1,[1,1])



C_FT = np.ones(len(E_C))
psi0 = np.ones(len(E_C))

#plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0,1,len(E_C))))
ax1 = plt.subplot(gs[1])

r_0 = 1e-20


#T = 0.2
#N = 50
#r = np.linspace(0,N*0.2,N,endpoint=False)

#k = fftfreq(N,T)[:N//2]

psi_FT = np.ones((len(k),len(E_C)))

for j in range(len(E_C)):
    
    print(j)
    
    [norm_int,res_int] = np.sqrt(scinteg.quad(BuschFunc.psi_s_int, 0,50, args=(E_C[j])))
    
    for i in range(len(k)):
        
        psi_FT[i,j] = scinteg.quad(amp_ft_integ, 0, 5, args=(k[i],E_C[j]))[0] 
        psi_FT[i,j] += scinteg.quad(amp_ft_integ, 5, 10, args=(k[i],E_C[j]))[0]
        psi_FT[i,j] = psi_FT[i,j]/norm_int
       # psi_FT[i,j] = 
        
       
    C_FT[j] = curve_fit(high_mom,k,np.abs(psi_FT[:,j])**2,p0=C_ad[j])[0]
    
    
    
  #  psi_FT[j] = fft(r**2*BuschFunc.psi_s(r,E_C[j])/norm_int)
    psi0[j] = (r_0*BuschFunc.psi_s(r_0,E_C[j])/norm_int)**2
    
  #  C_FT[j] = curve_fit(high_mom,k,np.abs(psi_FT[j])**2,p0=C_ad[j])[0]
    
    
    
 
    
    if np.any(np.abs(psi_FT[:,j])**2 < 1e-15):
        
        print(j)
  
ax1.plot(k,np.abs(psi_FT)**2)
#ax1.set_xticks([], minor=True)

#ax1.set_xticks([100, 150, 200, 250])
#ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

plt.xscale('log')
plt.yscale('log')

ax1.minorticks_off()
ax1.set_xticks([], minor=True)

plt.setp(ax1.spines.values(), linewidth=spine_width)
ax1.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)



#plt.xticks(fontsize=fontsize_ticks,**hfont)
plt.yticks(fontsize=fontsize_ticks,**hfont)

#plt.xlim([100,250])

#plt.xticks([100,150, 200, 250],fontsize=fontsize_ticks,**hfont)
#plt.yticks(fontsize=fontsize_ticks,**hfont)

plt.xlabel(r'$k$ $(a_{\mathrm{ho}}^{-1})$',fontsize=fontsize,**hfont)
plt.ylabel(r'$n(k)$ (arb.)',fontsize=fontsize,**hfont)

fig_mom.tight_layout(pad=1.5) 

#plt.savefig('Plots/Ch2_momentum_figure_plot.png')   
        
#plt.rcParams["axes.prop_cycle"] = plt.cycler('tab10')

fig = plt.figure(dpi=350,figsize=(6, 4))

plt.plot(E_C,C_ad,marker='o',linestyle='none')
plt.plot(E_C,psi0*16*np.pi**2,'.',linewidth=2)
plt.plot(E_C,C_FT,marker='o',linestyle='none')
#plt.ylim([0,8])

C_array = np.ones((len(E_C),5))
C_array[:,0] = E_C
C_array[:,1] = C_ad
C_array[:,2] = C_FT
C_array[:,3] = psi0



np.savetxt('C_array_lower.csv',C_array,delimiter=',')
        

        
        
    
    


