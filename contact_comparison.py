#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 00:27:19 2022

@author: robynlearn
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

import BuschFunc

#plot settings
hfont = {'fontname':'Helvetica'}
fontsize = 16
fontsize_ticks = 12
fontsize_leg = 12
spine_width = 1.5

hbar = 6.626e-34/(2*np.pi)
a_B = 5.29177e-11 #bohr radius [m]
m = 39.96399848*1.66053906660e-27

#original a0
C_array = np.loadtxt('C_array_orig_a0.csv',delimiter=",",dtype='float')
bound_array = np.loadtxt('bound_array_orig_a0.csv',delimiter=",",dtype='float')

plot_orig_a0=0
plot_B=1
plot_kF = 0

if plot_orig_a0 == 1:
    E_bnd_deep = -340
    
    fig = plt.figure(dpi=350,figsize=(8,4))
    #gs = gridspec.GridSpec(3,3)
    gs = gridspec.GridSpec(1,2)
    
    E_C = C_array[:,0]
    C_ad = C_array[:,1]
    C_FT = C_array[:,2]
    psi0 = C_array[:,3]
    
    E_ov = bound_array[:,0]
    overlap = bound_array[:,1]
    overlap_anal = bound_array[:,2]
    overlap_anal_deep = bound_array[:,3]
    
    ax0 = plt.subplot(gs[0,0])
    ax0.plot(E_C,C_FT,marker='o',linestyle='none')
    
    plt.xlim([1.5,3.5])
    plt.ylim([0,10])
    plt.xticks(fontsize=fontsize_ticks,**hfont)
    plt.yticks(fontsize=fontsize_ticks,**hfont)
    
    plt.setp(ax0.spines.values(), linewidth=spine_width)
    ax0.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)
    
    plt.xlabel(r'Energy $(\hbar \omega)}$',fontsize=fontsize,**hfont)
    plt.ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)',fontsize=fontsize,**hfont) 
    
    #ax1 = plt.subplot(gs[0])
    ax0.plot(E_C,C_ad,marker='o',linestyle='none')
    
    ax0.plot(E_C,psi0*8*np.sqrt(2)*np.pi**2,linewidth=2)
    ax0.plot(E_C,psi0*16*np.pi**2,'--',linewidth=2)
    
    plt.xlim([1.5,3.5])
    plt.ylim([0,10])
    plt.xticks(fontsize=fontsize_ticks,**hfont)
    plt.yticks(fontsize=fontsize_ticks,**hfont)
    
    
    ax2 = plt.subplot(gs[1])
    ax2.plot(E_ov,overlap**2,linewidth=2)
    ax2.plot(E_ov,overlap_anal**2,'.',markersize=6)
    #plt.plot(E_ov,overlap_anal**2*8*np.pi/BuschFunc.a_0_func(-200),linestyle='none',marker='o')
    
    plt.xlim([1.5,3.5])
    plt.ylim([0,0.05])
    plt.xticks(fontsize=fontsize_ticks,**hfont)
    plt.yticks(fontsize=fontsize_ticks,**hfont)
    
    plt.xlabel(r'Energy $(\hbar \omega)}$',fontsize=fontsize,**hfont)
    plt.ylabel(r'$N_{\mathrm{mol}}$ (arb.)',fontsize=fontsize,**hfont)
    
    plt.setp(ax2.spines.values(), linewidth=spine_width)
    ax2.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False) 
    
    fig.tight_layout(pad=1)
    #plt.savefig('Plots/Ch2_contact_orig_a0_plots.png')

#corrected a0
C_array = np.loadtxt('C_array.csv',delimiter=",",dtype='float')
bound_array = np.loadtxt('bound_array.csv',delimiter=",",dtype='float')

C_0_array = np.loadtxt('C_array_lower.csv',delimiter=",",dtype='float')



E_C = C_array[:,0]
C_ad = C_array[:,1]
C_FT = C_array[:,2]
psi0 = C_array[:,3]

E_C_0 = C_0_array[:,0]
C_ad_0 = C_0_array[:,1]
C_FT_0 = C_0_array[:,2]
psi0_0 = C_0_array[:,3]


if plot_B==1:

    V_L_list = np.array([50,100,200,300])
    
    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    hbar = 6.626*10**-34/(2*np.pi)
    m = 39.96399848*1.66054e-27
    E_R = hbar**2*k_L**2/(2*m)
    
    fig, axs = plt.subplots(ncols=2,dpi=600,figsize=(9,4))
    
    ax = axs[0]
    ax0 = axs[1]
    
    for V_L in V_L_list:
        
      #  V_L = V_L_list[i]
    #
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
        
        B = BuschFunc.B_func(BuschFunc.a_0_func(E_C), V_L)
        B_reff = BuschFunc.B_func_r_eff(E_C, V_L)
        
        B_reff_0 = BuschFunc.B_func_r_eff(E_C_0, V_L)
    
    
        ax.plot(B_reff,C_ad/a_ho*a_B,linestyle='none',marker='o',label=r'$V_{\mathrm{L}}$=%0.0f$E_{\mathrm{R}}$' % (V_L))
       # ax.plot(B_reff,BuschFunc.contact_upper_interp(B_reff, V_L)/a_ho*a_B)
        
        
    
        ax0.plot(B_reff_0,C_ad_0/a_ho*a_B,linestyle='none',marker='o',label=r'$V_{\mathrm{L}}$=%0.0f$E_{\mathrm{R}}$' % (V_L))
      #  ax0.plot(B_reff_0,BuschFunc.contact_lower_interp(B_reff_0, V_L)/a_ho*a_B)
      #  ax.plot(B_reff,C_ad/a_ho*a_B,linestyle='-',label=r'$V_{\mathrm{L}}$=%0.0f$E_{\mathrm{R}}$' % (V_L))
    
        
      #  ax.plot(B_reff,C_ad/a_ho,linestyle='-',label=r'$V_{\mathrm{L}}$=%0.0f$E_{\mathrm{R}}$' % (V_L))
    
    ax.set_xlim([185,250])
   # ax0.set_xlim([195,220])
    ax.set_ylim([-0.001,0.01]) 
  #  ax0.set_ylim([-0.001,0.01])
        
    
        
       # ax.plot(B_reff[np.argmax(C_ad)]*np.ones(10),np.arange(10)*10**7,'--',label='Max. @ %0.2f G' % (B_reff[np.argmax(C_ad)]))
       # ax.plot(B_reff[np.argmin(C_ad)]*np.ones(10),np.arange(10)*10**7,'--',label='Min. @ %0.2f G' % (B_reff[np.argmin(C_ad)]))
    
    ticks = np.arange(195,211,2.5)
    ticks0 = np.arange(200,211,2.5)
    ax.set_xticks(ticks)
 #   ax.set_yticks(ticks)
    
    ax.set_xticklabels(list(map(str,ticks)),fontsize=fontsize_ticks,**hfont)
   # ax.set_yticklabels(ax.get_yticklabels,fontsize=fontsize_ticks,**hfont)
    
    ax0.set_xticks(ticks0)
   # ax.set_yticks(ticks)
    
    ax0.set_xticklabels(list(map(str,ticks0)),fontsize=fontsize_ticks,**hfont)
 #   ax.set_yticklabels(fontsize=fontsize_ticks,**hfont)
    
    plt.setp(ax.spines.values(), linewidth=spine_width)
    ax.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)
    ax.grid('on')
    
    plt.setp(ax0.spines.values(), linewidth=spine_width)
    ax0.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)
    ax0.grid('on')
    
    ax.set_xlabel(r'B (G)',fontsize=fontsize,**hfont)
    #plt.ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)',fontsize=fontsize,**hfont) 
    ax.set_ylabel(r'Contact ($a_{\mathrm{B}}^{-1}$)',fontsize=fontsize,**hfont) 
    
    ax0.set_xlabel(r'B (G)',fontsize=fontsize,**hfont)
    #plt.ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)',fontsize=fontsize,**hfont) 
    ax0.set_ylabel(r'Contact ($a_{\mathrm{B}}^{-1}$)',fontsize=fontsize,**hfont) 
    plt.legend()
    
    ax.set_title('Upper branch')
    ax0.set_title('Lower branch')
    
    fig.tight_layout(pad=1)
    
    
    
    plt.savefig('Contact_01.png')




if plot_kF == 1:

    E_ov = bound_array[:,0]
    overlap = bound_array[:,1]
    overlap_anal = bound_array[:,2]
    overlap_anal_deep = bound_array[:,3]
    
    fig = plt.figure(dpi=600,figsize=(8,4))
    gs = gridspec.GridSpec(1,2)
    
    ax3 = plt.subplot(gs[0,0])
    ax3.plot(E_C,C_FT,marker='o',linestyle='none')
    ax3.plot(E_ov,overlap_anal**2*350,linestyle='none',marker='o')
    #ax3.plot(E_C,C_ad,marker='o')
    ax3.plot(E_C,C_ad,'o')
    ax3.plot(E_C,psi0*16*np.pi**2,'.',linewidth=2)
    ax3.plot(E_ov,overlap**2*260,'--',linewidth=3)
    
    #plt.xlim([1.5,3.5])
    #plt.ylim([0,10])
    plt.xticks(fontsize=fontsize_ticks,**hfont)
    plt.yticks(fontsize=fontsize_ticks,**hfont)
    
    plt.setp(ax3.spines.values(), linewidth=spine_width)
    ax3.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)
    
    plt.xlabel(r'Energy $(\hbar \omega)}$',fontsize=fontsize,**hfont)
    plt.ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)',fontsize=fontsize,**hfont) 
    
    
    
    ax4 = plt.subplot(gs[0,1])
    
    kF = np.sqrt(6)
    
    N = 1
    
    E_C =(BuschFunc.a_0_func(E_C))**(-1)/kF
    C_ad = C_array[:,1]/kF/N
    C_FT = C_array[:,2]/kF/N
    psi0 = C_array[:,3]/kF/N
    
    E_ov =(BuschFunc.a_0_func(E_ov))**(-1)/kF
    overlap = bound_array[:,1]/np.sqrt(kF*N)
    overlap_anal = bound_array[:,2]/np.sqrt(kF*N)
    overlap_anal_deep = bound_array[:,3]/np.sqrt(kF*N)
    
    ax4.plot(E_C,C_FT,marker='o',linestyle='none',label='Mom. tail')
    ax4.plot(E_ov,overlap_anal**2*350,linestyle='none',marker='o',label='Dim. overlap')
    #ax3.plot(E_C,C_ad,marker='o')
    ax4.plot(E_C,C_ad,'o',label='Adiabatic')
    ax4.plot(E_C,psi0*16*np.pi**2,'.',linewidth=2,label='Short-range')
    ax4.plot(E_ov,overlap**2*260,'--',linewidth=3,label='Int. overlap')
    
    plt.xlim([-5,5])
    #plt.ylim([2.8,3.]/2)
    plt.xticks(fontsize=fontsize_ticks,**hfont)
    plt.yticks(fontsize=fontsize_ticks,**hfont)
    
    plt.setp(ax4.spines.values(), linewidth=spine_width)
    ax4.tick_params(which='both',direction='in', length=5, width=spine_width,top=False,right=False)
    
    plt.xlabel(r'Inverse scattering length $(k_{\mathrm{F}})}$',fontsize=fontsize,**hfont)
    plt.ylabel(r'Contact ($k_{\mathrm{F}}$)',fontsize=fontsize,**hfont) 
    plt.legend()
    
    fig.tight_layout(pad=1)
    #plt.savefig('Plots/Ch2_contact_summary_plots.png')


