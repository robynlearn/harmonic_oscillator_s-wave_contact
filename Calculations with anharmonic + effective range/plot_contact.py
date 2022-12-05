#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 12:58:45 2022

@author: robynlearn
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.io as scio

import BuschFunc


###### constants #######

hbar = 6.626e-34/(2*np.pi)
a_B = 5.29177e-11 #bohr radius [m]
m = 39.96399848*1.66053906660e-27 # [kg]

r_eff = 98*a_B

lambda_L = 1054e-9 #laser wavelength [m]
k_L = 2*np.pi/lambda_L
E_R = hbar**2*k_L**2/(2*m)

B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
a_bg = 166.978*a_B #background scattering rate
del_B = 6.910 #resonance width [G]

mu_B = 9.2740100783e-24/10000 #J/G
mu_mol = 1.68*mu_B

###### flags #######

do_gamma = 1
do_gammakHz = 0
plot_lifetime_powers = 0
plot_CB = 0
plot_frac_C = 0
plot_EB = 0




##### Do the stuff #######


if do_gamma == 1:
    
    interp = 1
    
    #guesstimate molecular rabi freq plus lifetime
    
  #  gamma_mol = 2*np.pi*10*1e6
  #  rabi_mol = 2*np.pi*1*1e6
    
    gamma_mol = 2*np.pi*26*1e6
    rabi_mol = 2*np.pi*2.6*1e6
    
    data = scio.loadmat('/Volumes/GoogleDrive/My Drive/Lattice Shared/SharedData/2022 PA experiment/11_25 lattice_lifetime/lifetime.mat')
    data = data['lifetime']
    B_data = data[0][0][0][:,0] #G
    lifetimes = (data[0][0][1][0])*1000 #1/s
    lifetimes_err = (data[0][0][2][0])*1000 #1/s
    
    B = np.linspace(200,210,200)
    
    #V_L_array = np.array([50,100,200,300])
    V_L_array = [200]#[200]
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for V_L in V_L_array:
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
           
      #  a = BuschFunc.a_interp_B_97_1st_branch(B, BuschFunc.E_interp_1st_branch_97(B, V_L), V_L)*a_ho
       # dadB = np.gradient(a**-1,B)/a_ho
       
        if interp == 1:
            
            C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho#*a_B
            C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho#*a_B
        
        else:
            
            C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
            B = BuschFunc.B_func_97(C_array[:,0], V_L)
            C_lower = C_array[:,2]/a_ho
        
        ax.errorbar(B_data,lifetimes,yerr=lifetimes_err,marker='o',linestyle='none')
        
        R = hbar**2/(m*a_bg*mu_mol*del_B)
        
        N_mol_upper = C_upper*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)
       # N_mol = C_lower*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)/2
        N_mol = C_lower*del_B**2*R/(4*np.pi*((B - B_0 - del_B )**2+(rabi_mol**2/(2*gamma_mol*del_B*mu_mol/hbar))**2))/2
     #   N_mol = C_lower*R/(4*np.pi)*(1-a_bg/a)**2
       # N_mol = -C_lower*hbar**2/(mu_mol*4*np.pi*m)*dadB
        
        gamma_upper = 2*N_mol_upper*rabi_mol**2/(gamma_mol)
        gamma = 2*N_mol*rabi_mol**2/(gamma_mol)
            
       # ax.plot(B,frac,label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)
       # ax.plot(B,C_upper)
       # ax.plot(B,C_lower)
        ax.plot(B,gamma,'.')
        #ax.plot(B,gamma_upper,'.')
        
      #  ax.plot(B_res*np.ones(10),np.linspace(-0.1,1.3,10),'--',color='gray')
        
      #  plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
        
    #ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([200,210])
    plt.ylim([0,0.1e6])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'$\Gamma$ (rad)')
    
    plt.legend()
    
if do_gammakHz == 1:
    
    interp = 1
    
    #guesstimate molecular rabi freq plus lifetime
    
  #  gamma_mol = 2*np.pi*10*1e6
  #  rabi_mol = 2*np.pi*1*1e6
    
    gamma_mol = 2*np.pi*26*1e6 #max power 
    rabi_mol = 2*np.pi*2.6*1e6 #max power
    
    data = scio.loadmat('/Volumes/GoogleDrive/My Drive/Lattice Shared/SharedData/2022 PA experiment/11_25 lattice_lifetime/lifetime.mat')
    data = data['lifetime']
    B_data = data[0][0][0][:,0] #G
    gamma_data = (data[0][0][1][0]) #1/ms
    gamma_err = (data[0][0][2][0]) #1/ms
    
    lifetimes = 1/gamma_data #ms
    lifetimes_err = gamma_err/gamma_data**2 #ms
    
    
    B = np.linspace(200,210,200)
    
    #V_L_array = np.array([50,100,200,300])
    V_L_array = [200]
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for V_L in V_L_array:
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
           
      #  a = BuschFunc.a_interp_B_97_1st_branch(B, BuschFunc.E_interp_1st_branch_97(B, V_L), V_L)*a_ho
       # dadB = np.gradient(a**-1,B)/a_ho
       
        if interp == 1:
            
            C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho#*a_B
            C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho#*a_B
        
        else:
            
            C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
            B = BuschFunc.B_func_97(C_array[:,0], V_L)
            C_lower = C_array[:,2]/a_ho
        
        ax.errorbar(B_data,lifetimes,yerr=lifetimes_err,marker='o',linestyle='none')
        
        R = hbar**2/(m*a_bg*mu_mol*del_B)
        
        N_mol_upper = C_upper*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)
       # N_mol = C_lower*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)/2
        N_mol = C_lower*del_B**2*R/(4*np.pi*((B - B_0 - del_B )**2+(rabi_mol**2/(2*gamma_mol*del_B*mu_mol/hbar))**2))/2
     #   N_mol = C_lower*R/(4*np.pi)*(1-a_bg/a)**2
       # N_mol = -C_lower*hbar**2/(mu_mol*4*np.pi*m)*dadB
        
        gamma_upper = 2*N_mol_upper*rabi_mol**2/(gamma_mol)
        gamma = 2*N_mol*rabi_mol**2/(gamma_mol)
            
       # ax.plot(B,frac,label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)
       # ax.plot(B,C_upper)
       # ax.plot(B,C_lower)
        ax.plot(B,1/gamma*1000,'.')
        #ax.plot(B,gamma_upper,'.')
        
      #  ax.plot(B_res*np.ones(10),np.linspace(-0.1,1.3,10),'--',color='gray')
        
      #  plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
        
    #ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([200,210])
   # plt.ylim([0,0.1e6])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'$\tau$ (ms)')
    
    plt.legend()
    
if plot_lifetime_powers == 1:
    
    interp = 1
    
    #guesstimate molecular rabi freq plus lifetime
    
  #  gamma_mol = 2*np.pi*10*1e6
  #  rabi_mol = 2*np.pi*1*1e6
    
    gamma_mol = 2*np.pi*26*1e6 #max power
    rabi_mol_max = 2*np.pi*2.6*1e6 #max power
    
    P_max = 1.15
    P_array = np.array([0.03,0.2,0.6,1.15])
    
    data = scio.loadmat('/Volumes/GoogleDrive/My Drive/Lattice Shared/SharedData/2022 PA experiment/11_25 lattice_lifetime/lifetime.mat')
    data = data['lifetime']
    B_data = data[0][0][0][:,0] #G
    gamma_data = (data[0][0][1][0]) #1/ms
    gamma_err = (data[0][0][2][0]) #1/ms
    
    lifetimes = 1/gamma_data #ms
    lifetimes_err = gamma_err/gamma_data**2 #ms
    
    
    B = np.linspace(200,210,200)
    
    #V_L_array = np.array([50,100,200,300])
    V_L = 200#np.array([200])
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for P in P_array:
        
        rabi_mol = rabi_mol_max*np.sqrt(P/P_max)
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
           
      #  a = BuschFunc.a_interp_B_97_1st_branch(B, BuschFunc.E_interp_1st_branch_97(B, V_L), V_L)*a_ho
       # dadB = np.gradient(a**-1,B)/a_ho
       
        if interp == 1:
            
            C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho#*a_B
            C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho#*a_B
        
        else:
            
            C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
            B = BuschFunc.B_func_97(C_array[:,0], V_L)
            C_lower = C_array[:,2]/a_ho
        
       # ax.errorbar(B_data,lifetimes,yerr=lifetimes_err,marker='o',linestyle='none')
        
        R = hbar**2/(m*a_bg*mu_mol*del_B)
        
        N_mol_upper = C_upper*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)
       # N_mol = C_lower*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)/2
        N_mol = C_lower*del_B**2*R/(4*np.pi*((B - B_0 - del_B )**2+(rabi_mol**2/(2*gamma_mol*del_B*mu_mol/hbar))**2))/2
     #   N_mol = C_lower*R/(4*np.pi)*(1-a_bg/a)**2
       # N_mol = -C_lower*hbar**2/(mu_mol*4*np.pi*m)*dadB
        
        gamma_upper = 2*N_mol_upper*rabi_mol**2/(gamma_mol)
        gamma = 2*N_mol*rabi_mol**2/(gamma_mol)
            
       # ax.plot(B,frac,label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)
       # ax.plot(B,C_upper)
       # ax.plot(B,C_lower)
        ax.plot(B,1/gamma*1000,'-',label='% .2f mW' % P)
        #ax.plot(B,gamma_upper,'.')
        
      #  ax.plot(B_res*np.ones(10),np.linspace(-0.1,1.3,10),'--',color='gray')
        
      #  plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
        
    ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([202,209])
    plt.ylim([0,50])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'$\tau$ (ms)')
    
    plt.legend()
    
if plot_CB == 1:
    
    interp = 1
    
   # B = np.linspace(200,210,1000)
    B = np.linspace(200,210,1000)
    
    
   # V_L_array = np.array([50,100,200,300])
    V_L_array = [200]
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for V_L in V_L_array:
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
        
        if interp == 1:
            
            B_lower = B
            B_upper = B
            
            C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho#*a_B
            C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho#*a_B

        else:
            
            upper_B = 260
            lower_B = 200
            
            C_array_upper = np.loadtxt('C_array_upper_'+str(V_L)+'ER.csv',delimiter=',')
            B_full_upper = BuschFunc.B_func_97(C_array_upper[:,0], V_L)
            #B_upper = B_full_upper
            B_upper = B_full_upper[(B_full_upper>lower_B) & (B_full_upper<upper_B)]
            C_upper = C_array_upper[:,2]/a_ho
            C_upper = C_upper[(B_full_upper>lower_B) & (B_full_upper<upper_B)]
            
            C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
            B_full = BuschFunc.B_func_97(C_array[:,0], V_L)
            B_lower = B_full[(B_full>lower_B) & (B_full<upper_B)]
            C_lower = C_array[:,2]/a_ho
            C_lower = C_lower[(B_full>lower_B) & (B_full<upper_B)]
            

        
        ax.plot(B_upper,C_upper,'.',label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)
        ax.plot(B_lower,C_lower,'.',label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)
        
        
      #  plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
        
    #ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([200,210])
  #  plt.ylim([0,1e9])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'$C$ ($m^{-1}$)')
    
    plt.legend()
    
    
if plot_frac_C == 1:
    
    B = np.linspace(200,210,1000)
    
    V_L_array = [200]
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for V_L in V_L_array:
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
           
        C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho
        C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho
        
        frac = C_upper/C_lower
        
            
        ax.plot(B,frac,'.',label='$V_{\mathrm{L}}$=%.0f $E_{\mathrm{R}}$' % V_L)

        
    #ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([200,210])
    plt.ylim([-0.001,2])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'C($\beta=1$)/C($\beta=0$)')
    
    plt.legend()
    
if plot_EB == 1:
    
   # B = np.linspace(200,210,200)
    
    
    V_L_array = np.array([50,100,200,300])
   # V_L_array = [200]
    
    fig, ax = plt.subplots(dpi=350,figsize=(4,3))
    
    for V_L in V_L_array:
        
        omega = 2*E_R*np.sqrt(V_L)/hbar
        a_ho = np.sqrt(2*hbar/(m*omega))
    
       # C_upper = BuschFunc.contact_upper_interp(B,V_L)/a_ho*a_B
       # C_lower = BuschFunc.contact_lower_interp(B,V_L)/a_ho#*a_B
       
        C_array_upper = np.loadtxt('C_array_upper_'+str(V_L)+'ER.csv',delimiter=',')
        B_full_upper = BuschFunc.B_func_97(C_array_upper[:,0], V_L)
        E_C_upper = C_array_upper[:,1]
        C_upper = C_array_upper[:,2]/a_ho
        
        B_upper = B_full_upper
      #  B_upper = B_full_upper[(B_full_upper>200) & (B_full_upper<210)]
      #  E_C_upper = E_C_upper[(B_full_upper>200) & (B_full_upper<210)]
      #  C_upper = C_upper[(B_full_upper>200) & (B_full_upper<210)]
        
        
        
        C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
        B_full = BuschFunc.B_func_97(C_array[:,0], V_L)
        E_C = C_array[:,1]
        C_lower = C_array[:,2]/a_ho
        
        B = B_full
       # B = B_full[(B_full>200) & (B_full<210)]
       # E_C = E_C[(B_full>200) & (B_full<210)]
       # C_lower = C_lower[(B_full>200) & (B_full<210)]
        
        ax.plot(B,E_C,'.')
        ax.plot(B_upper,E_C_upper,'.')
        
        
        
      #  plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
        
    #ax.set_yscale('log')   
        
    ax.grid('on')
    plt.xlim([195,260])
    plt.ylim([-3,4])

    plt.title('$V_{\mathrm{L}}$=%.0f' % V_L)
    
    
    plt.xlabel('B field (G)')
    plt.ylabel(r'Energy ($\hbar \omega$)')
    
    plt.legend()