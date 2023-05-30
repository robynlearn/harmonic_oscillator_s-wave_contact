#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:40:08 2023

@author: robynlearn
"""

import numpy as np
import BuschFunc
import matplotlib.pyplot as plt

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


B = 196
V_L = 100

shift = BuschFunc.E_interp_1st_branch_95(B,V_L) - BuschFunc.E_interp_2nd_branch_97(B,V_L)

omega = 2*E_R*np.sqrt(V_L)/hbar

shift_kHz = shift*omega/(1000*2*np.pi)

print(shift_kHz)

B = np.array([204])
shiftcalc = BuschFunc.rf_spec_shift(B, V_L)


print(shiftcalc)

B = np.linspace(195,210,100)

fig = plt.figure(figsize=[4,3],dpi=350)
plt.plot(B,BuschFunc.rf_spec_shift(B, V_L),'.')

plt.xlabel('B field (G)')
plt.ylabel('75 shift (kHz)')