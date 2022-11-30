#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:09:10 2022

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


######### Define functions ###########

#scattering length for given energy
def a_0_func(E):
    
  #  a_0 = 1/np.sqrt(2)*scsp.gamma(-E/2+1/4)/scsp.gamma(-E/2+3/4)
   # a_0 = scsp.gamma(-E/2+1/4)/scsp.gamma(-E/2+3/4)
   
    a_0 = 1/2*scsp.gamma(-E/2+1/4)/scsp.gamma(-E/2+3/4)
    
    return a_0

#B-field dependence of scattering length for 40K -7,-9 Feschbach resonance

def B_func(a_0,V_L):
    
    B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
    a_B = 5.29177e-11 #bohr radius [m]
    a_bg = 166.978*a_B #background scattering rate
    del_B = 6.910 #resonance width [G]

  #  V_L = 200 #lattice depth [E_R]
    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L

    a = a_0/(k_L*V_L**(1/4)) #scattering length in metres for given lattice depth


    B = del_B/(1 - a/a_bg) + B_0
    
    return B

def B_func_r_eff(E,V_L):
    
    a_0 = a_0_func(E)
    
    hbar = 6.626e-34/(2*np.pi)
    a_B = 5.29177e-11 #bohr radius [m]
    m = 39.96399848*1.66053906660e-27
    
    
    B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
    a_bg = 166.978*a_B #background scattering rate
    del_B = 6.910 #resonance width [G]

  #  V_L = 200 #lattice depth [E_R]
    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    
    
    hbar_omega = 2*np.sqrt(V_L)*E_R

    a = a_0/(k_L*V_L**(1/4)) #scattering length in metres for given lattice depth
    
    r_eff = 98*a_B
    B_0 = B_0 + del_B/(1 - hbar**2/(r_eff*E*hbar_omega*(m/2)*a_bg))
    
    B = del_B/(1 - a/a_bg) + B_0
    
    return B
    

def B_func_95(a_0,V_L):
    
    B_0 = 224.2 #95 Feshbach resonance [G] (chip lab PRR values)
    a_B = 5.29177e-11 #bohr radius [m]
    a_bg = 167.3*a_B #background scattering rate
    del_B = 7.2 #resonance width [G]

   # V_L = 200 #lattice depth [E_R]
    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L

    a = a_0/(k_L*V_L**(1/4)) #scattering length in metres for given lattice depth

    B = del_B/(1 - a/a_bg) + B_0
    
    return B

def E_interp_1st_branch_95(B):
    
    E = np.linspace(1.5+0.00000001,2.5-0.00000001,100)
#     E = np.linspace(-50,1.5-0.000001,100)
    a_0 = a_0_func(E)
    B_interp = B_func_95(a_0)
    
    E_interp_func = scinterp.interp1d(B_interp,E)
    E_interp = E_interp_func(B)
    
    return E_interp

def E_interp_1st_branch_97(B,V_L):
    
    E = np.linspace(-10,1.5+0.14,200)

    B_interp = B_func_r_eff(E,V_L)
    
   
    E_interp_func = scinterp.interp1d(B_interp,E)
    E_interp = E_interp_func(B)
    
    return E_interp

def E_interp_2nd_branch_97(B,V_L):
    
    E = np.linspace(1.8,3.5-0.00000001,100)

    B_interp = B_func_r_eff(E,V_L)
    
   
    E_interp_func = scinterp.interp1d(B_interp,E)
    E_interp = E_interp_func(B)
    
    return E_interp

def E_interp_a(a):
    
    E = np.linspace(1.5+0.00000001,3.5-0.00000001,100)
    
    a_0 = a_0_func(E)
    
    E_interp_func = scinterp.interp1d(a_0,E)
    
    E_interp = E_interp_func(a)
    
    return E_interp


#harmonic oscillator energy for given n
def E_n(n):  
    
    E_n = 2*n + 3/2
    
    return E_n

#creates array of harm. osc. energy by column, for plotting against a_0
def E_n_array(n,dim):
   
   E_n_array = np.ones([dim,n])
   
   for i in range(n):
       E_n_array[:,i] = E_n(i)
       
   return E_n_array
   
#wavefunction of interacting state (not normalized)

def psi_s(r,E):
    
    nu = E/2 - 3/4
    
    psi_s = 0.5*np.pi**(-3/2)*np.exp(-r**2/2)*scsp.gamma(-nu)*scsp.hyperu(-nu,3/2,r**2)
    
    return psi_s

#harmonic oscillator wavefunction
def psi_harm(r,n):
    
    psi = np.pi**(-3/4)*(scsp.genlaguerre(n,0.5)(0))**(-1/2)*np.exp(-r**2/2)*(scsp.genlaguerre(n,0.5)(r**2))

    return psi

#bound state wavfunction (not normalized)

def psi_bnd(r,E):
    
    psi = 1/r*np.exp(-r/a_0_func(E))
    
    return psi

def psi_harm_p(r,n):
    
    psi = r*np.exp(-r**2/2)*(scsp.genlaguerre(n,1.5)(r**2))
    
    return psi


#interacting state amplitude integrand

def psi_s_int(r,E):
    
    nu = E/2 - 3/4
    
    psi_amp = 4*np.pi*r**2*(np.abs(0.5*np.pi**(-3/2)*np.exp(-r**2/2)*scsp.gamma(-nu)*scsp.hyperu(-nu,3/2,r**2)))**2
    
    return psi_amp

def psi_bnd_anal_int(r,E):
    
    psi_amp = 4*np.pi*r**2*np.abs(psi_bnd(r,E))**2
    
    return psi_amp

    #harm osc amplitude integrand

def psi_harm_int(r,n):
    
    psi_amp = 4*np.pi*r**2*(np.abs(np.pi**(-3/4)*(scsp.genlaguerre(n,0.5)(0))**(-1/2)*np.exp(-r**2/2)*(scsp.genlaguerre(n,0.5)(r**2))))**2

    return psi_amp 

def psi_harm_p_int(r,n):
    
    psi_amp = 4*np.pi*r**2*(np.abs(r*np.exp(-r**2/2)*(scsp.genlaguerre(n,1.5)(r**2))))**2
    
    return psi_amp

#overlap integrals

def overlap_int(r,n,E):
    
    nu = E/2 - 3/4
    
    overlap = 4*np.pi*r**2*0.5*np.pi**(-3/2)*np.exp(-r**2/2)*scsp.gamma(-nu)*scsp.hyperu(-nu,3/2,r**2)*np.pi**(-3/4) \
    *(scsp.genlaguerre(n,0.5)(0))**(-1/2)*np.exp(-r**2/2)*(scsp.genlaguerre(n,0.5)(r**2))
    
    return overlap

def overlap_int_bnd(r,E,E_bnd):
    
    overlap = 4*np.pi*r**2*psi_s(r,E)*psi_s(r,E_bnd)
    
    return overlap

def overlap_int_bnd_anal(r,E,E_bnd):
    
    overlap = 4*np.pi*r**2*psi_s(r,E)*psi_bnd(r,E_bnd)
    
    return overlap

def contact_upper_interp(B,V_L):
    
    C_array = np.loadtxt('C_array.csv',delimiter=",",dtype='float')
    
    E_C = C_array[:,0]
    C_ad = C_array[:,1]
   # C_FT = C_array[:,2]
   # psi0 = C_array[:,3]
    
    B_reff = B_func_r_eff(E_C, V_L)
    
    C_interp_func = scinterp.interp1d(B_reff,C_ad,'cubic')
    C_interp = C_interp_func(B)
    
    return C_interp

def contact_lower_interp(B,V_L):
    
    C_array = np.loadtxt('C_array_lower.csv',delimiter=",",dtype='float')
    
    E_C = C_array[:,0]
    C_ad = C_array[:,1]
  #  C_FT = C_array[:,2]
  #  psi0 = C_array[:,3]
    
    B_reff = B_func_r_eff(E_C, V_L)
    
    C_interp_func = scinterp.interp1d(B_reff,C_ad,'cubic')
    C_interp = C_interp_func(B)
    
    return C_interp
    
    
    
    
    
    












