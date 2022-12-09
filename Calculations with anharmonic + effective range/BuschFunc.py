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

#scattering length for given energy w. effective range
def a_0_func(E,V_L):
    
    a_B = 5.29177210903e-11 #bohr radius [m]
    hbar = 6.62607015e-34/(2*np.pi)
    m = 39.96399848*1.66053906660e-27

    lambda_L = 1054e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    
    omega = 2*E_R*np.sqrt(V_L)/hbar # [Hz]
    a_ho = np.sqrt(2*hbar/(m*omega)) #[m]
    
    r_eff = scsp.gamma(1/4)**2/(3*np.pi)*2*65.0223142660582*a_B #[m]
  #  r_eff = 98*a_B
    
    shift = E*r_eff/a_ho
   
    a_0 = 1/(2*scsp.gamma(-E/2+3/4)/scsp.gamma(-E/2+1/4) + shift)
    
   # a_0 = a_E_func(E)
    
    return a_0

#energy-dependent scattering length for given energy
def a_E_func(E):
    
    a_0 = 1/2*scsp.gamma(-E/2+1/4)/scsp.gamma(-E/2+3/4)
    
    return a_0

#B-field dependence of scattering length for 40K -7,-9 Feschbach resonance

def B_func_97(a_0,V_L):

    hbar = 6.62607015e-34/(2*np.pi)
    m = 39.96399848*1.66053906660e-27    
    a_B = 5.29177210903e-11 #bohr radius [m]
    
    B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
    a_bg = 166.978*a_B #background scattering rate
    del_B = 6.910 #resonance width [G]

    lambda_L = 1054e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    
    omega = 2*E_R*np.sqrt(V_L)/hbar    
    a_ho = np.sqrt(2*hbar/(m*omega))
  
    a = a_0*a_ho #[m]

    B = del_B/(1 - a/a_bg) + B_0
    
    return B
    

def B_func_95(a_0,V_L):
    
    hbar = 6.62607015e-34/(2*np.pi)
    a_B = 5.29177210903e-11 #bohr radius [m]
    m = 39.96399848*1.66053906660e-27 
    
    B_0 = 224.2 #95 Feshbach resonance [G] (chip lab PRR values)
    a_bg = 167.3*a_B #background scattering rate
    del_B = 7.2 #resonance width [G]


    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    
    omega = 2*np.sqrt(V_L)*E_R/hbar
    a_ho = np.sqrt(2*hbar**2/(m*hbar*omega))
    
    a = a_0*a_ho #[m]
    
    B = del_B/(1 - a/a_bg) + B_0

    return B


def E_interp_1st_branch_95(B,V_L):
    
    E = np.linspace(1.6,3.6,500)

    B_interp = B_func_95(E,V_L)
    
    E_interp_func = scinterp.interp1d(B_interp,E,kind='cubic')
    E_interp = E_interp_func(B)
    
    return E_interp

def E_interp_1st_branch_97(B,V_L):
    
    E = np.linspace(-10,1.7,500)

    B_interp = B_func_97(E,V_L)
   
    E_interp_func = scinterp.interp1d(B_interp,E,kind='cubic')
    E_interp = E_interp_func(B)
    
    return E_interp


def E_interp_2nd_branch_97(B,V_L):
    
    E = np.linspace(1.8,3.5-0.00000001,100)

    B_interp = B_func_97(E,V_L)
   
    E_interp_func = scinterp.interp1d(B_interp,E,kind='cubic')
    E_interp = E_interp_func(B)
    
    return E_interp

def E_interp_a(a):
    
    E = np.linspace(1.5+0.00000001,3.5-0.00000001,100)
    
    a_0 = a_0_func(E)
    
    E_interp_func = scinterp.interp1d(a_0,E,kind='cubic')
    
    E_interp = E_interp_func(a)
    
    return E_interp

def a_B_97_1st_branch(B,E,V_L):
    
    hbar = 6.62607015e-34/(2*np.pi)
    a_B = 5.29177210903e-11 #bohr radius [m]
    m = 39.96399848*1.66053906660e-27 
    
    B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
    a_bg = 166.978*a_B #background scattering rate
    del_B = 6.910 #resonance width [G]

    lambda_L = 1064e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    
    
    hbar_omega = 2*np.sqrt(V_L)*E_R
    a_ho = np.sqrt(2*hbar**2/(m*hbar_omega))
    
    a = a_bg*(1 - del_B/(B-B_0))/a_ho #[a_ho]
    
    return a

def a_interp_B_97_1st_branch(B,E,V_L):
    
    E_interp = np.linspace(-10,1.7,200)
    
    a_interp = a_0_func(E_interp)
    
    B_interp = B_func_97(E_interp,V_L)
    a_interp_func = scinterp.interp1d(B_interp,a_interp,kind='cubic')
    
    a_0 = a_interp_func(B)
    
    return a_0


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


##### various integrands for anharmonic shift ####

def psi_r4_int(r,E):
    
    integ = r**4*(psi_s(r,E))**2*4*np.pi
   # integ = r**4*(psi_harm(r,0))**2*4*np.pi
    
    return integ

def psi_r6_int(r,E):
    
    integ = r**6*(psi_s(r,E))**2*4*np.pi
   # integ = r**6*(psi_harm(r,0))**2*4*np.pi
    
    return integ



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

############ Overlap integrals ########

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

######## Contact functions #########

def calc_contact_lower(V_L):
    
    E_C = np.concatenate((np.arange(-3,0.6,0.001),np.arange(0.55,1.8,0.0002)))
    
    C_ad = -8*np.pi*np.gradient(E_C,1/a_0_func(E_C,V_L))
    
    return C_ad
    
    
    

def contact_upper_interp(B,V_L):
    
    C_array = np.loadtxt('C_array_upper_'+str(V_L)+'ER.csv',delimiter=',')
    
    a_C = C_array[:,0]
    C_ad = C_array[:,2]

    B_reff = B_func_97(a_C, V_L)
    
    C_interp_func = scinterp.interp1d(B_reff,C_ad,'cubic')
    C_interp = C_interp_func(B)
    
    return C_interp

def contact_lower_interp(B,V_L):
    
    C_array = np.loadtxt('C_array_lower_'+str(V_L)+'ER.csv',delimiter=',')
   
    a_C = C_array[:,0]
    C_ad = C_array[:,2]

    B_reff = B_func_97(a_C, V_L)

    C_interp_func = scinterp.interp1d(B_reff,C_ad,'cubic')
    C_interp = C_interp_func(B)
    
    return C_interp


##### anharmonic energy shift ####

def anharm_shift(E):
    
    [norm_int,res_int] = np.sqrt(scinteg.quad(psi_s_int, 0,20, args=(E)))
  #  [norm_int,res_int] = np.sqrt(scinteg.quad(psi_harm, 0,20, args=(0)))

    shift_4 = (1/8)*scinteg.quad(psi_r4_int,0,20,args=(E))[0]/norm_int**2
    shift_6 = (1/40)*scinteg.quad(psi_r6_int,0,20,args=(E))[0]/norm_int**2
    shift = (shift_6 + shift_4 + 3/32)
    

    
    
    return shift

def anharm_shift_E(E,V_L):
    
    hbar = 6.62607015e-34/(2*np.pi)
    m = 39.96399848*1.66053906660e-27

    lambda_L = 1054e-9 #laser wavelength [m]
    k_L = 2*np.pi/lambda_L
    E_R = hbar**2*k_L**2/(2*m)
    a_L = lambda_L/2
    
    omega = 2*E_R*np.sqrt(V_L)/hbar # [Hz]
    a_ho = np.sqrt(2*hbar/(m*omega)) #[m]
    
    [norm_int,res_int] = np.sqrt(scinteg.quad(psi_s_int, 0,10, args=(E)))

    shift_4 = scinteg.quad(psi_r4_int,0,100,args=(E))[0]/norm_int**2
    shift_6 = scinteg.quad(psi_r6_int,0,100,args=(E))[0]/norm_int**2
   # shift_E = -1*V_L*(E_R/(hbar*omega))*(a_ho*np.pi/a_L)**4*(shift_6 + shift_4 + 3/32)
    shift_E = -np.sqrt(V_L)/2*(a_ho*np.pi/a_L)**4*((1/40)*shift_6 + (1/8)*shift_4 + 3/32)
    
    
    return shift_E
    
    
    
    
    
    












