#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:24:40 2023

@author: robynlearn
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

import os, os.path

import scipy.special as scsp
import scipy.integrate as scinteg
import scipy.interpolate as scinterp
from scipy.optimize import curve_fit

import scipy.io as scio




import BuschFunc

hbar = 6.626e-34/(2*np.pi)
a_B = 5.29177e-11 #bohr radius [m]
m = 39.96399848*1.66053906660e-27


B_0 = 202.15 #Feshbach resonance [G] (chip lab PRR values)
a_bg = 166.978*a_B #background scattering rate
del_B = 6.910 #resonance width [G]

mu_B = 9.2740100783e-24/10000 #J/G
mu_mol = 1.68*mu_B
R = hbar**2/(m*a_bg*mu_mol*del_B)

#V_L = 200 #lattice depth [E_R]
lambda_L = 1064e-9 #laser wavelength [m]
k_L = 2*np.pi/lambda_L
E_R = hbar**2*k_L**2/(2*m)

def N_mol_func(B,C):
    
    N_mol = C*del_B**2*R/(4*np.pi*(B - B_0 - del_B )**2)
    
    return N_mol
    
    

#plot settings
hfont = {'fontname':'Helvetica'}
fontsize = 16
fontsize_ticks = 12
fontsize_leg = 12
spine_width = 1.5
linewidth=2



dim = 500
n = 5



##### Settings for Jose's calculations #####

freq = 128.66
#freq = 89.94
#freq = 69.93

omega = 2*np.pi*freq*10**3
V_L = (hbar*omega/(2*E_R))**2
a_ho = np.sqrt(2*hbar/(m*omega))
           

branch = 3 #3 is lower branch, 4 is upper branch
fort_num = 2000 + branch
plot_B = [190,202.075,209.075]

branch_upper = 4
fort_num_upper = 2000 + branch_upper


    
B = np.linspace(195.25,208.55,77)

cc_norm = np.ones(len(B))
oc_norm = np.ones(len(B))

oc_norm_shortrange = np.ones(len(B))
cc_norm_shortrange = np.ones(len(B))

cc_norm_turningpoint = np.ones(len(B))
oc_norm_turningpoint = np.ones(len(B))

cc_norm_upper = np.ones(len(B))
oc_norm_upper = np.ones(len(B))

V_L = 212.8074


C_array = np.loadtxt('/Users/robynlearn/Documents/GitHub/harmonic_oscillator_s-wave_contact/Calculations with anharmonic + effective range/C_array_lower_212.8074ER_no_anharm.csv',delimiter=',')
C_array_upper = np.loadtxt('/Users/robynlearn/Documents/GitHub/harmonic_oscillator_s-wave_contact/Calculations with anharmonic + effective range/C_array_upper_212.8074ER_no_anharm.csv',delimiter=',')


C_array = np.delete(C_array,slice(58024,58031),0)

B_upper = C_array_upper[:,3]
C_upper = C_array_upper[:,2]/a_ho

B_lower = C_array[:,3]
C_lower = C_array[:,2]/a_ho

C_array_no_eff_r = np.loadtxt('/Users/robynlearn/Documents/GitHub/harmonic_oscillator_s-wave_contact/Calculations with anharmonic + effective range/C_array_lower_212.8074ER_no_anharm_no_eff_r.csv',delimiter=',')
C_array_no_eff_r_upper = np.loadtxt('/Users/robynlearn/Documents/GitHub/harmonic_oscillator_s-wave_contact/Calculations with anharmonic + effective range/C_array_upper_212.8074ER_no_anharm_no_eff_r.csv',delimiter=',')

C_lower_no_eff_r = C_array_no_eff_r[:,2]/a_ho
B_lower_no_eff_r = C_array_no_eff_r[:,3]

C_upper_no_eff_r = C_array_no_eff_r_upper[:,2]/a_ho
B_upper_no_eff_r = C_array_no_eff_r_upper[:,3]




N_mol = N_mol_func(B_lower,C_lower)
N_mol_upper= N_mol_func(B_upper,C_upper)

N_mol_no_eff_r = N_mol_func(B_lower_no_eff_r,C_lower_no_eff_r)
N_mol_no_eff_r_upper = N_mol_func(B_upper_no_eff_r,C_upper_no_eff_r)

E_mol_lower = np.ones(len(B))
E_mol_upper = np.ones(len(B))
a = np.ones(len(B))

    
for i in range(len(B)):
    
    E_mol = np.loadtxt("/Users/UniversityThings/Thywissen Lab/MSc Project/Data-40K-ab/who="+("%.2f" % freq)+"kHz/Bfield="+("%.7f" % B[i])+"/Energies.data")
   # B = E[0]
    E_mol_lower[i] = E_mol[3]*4.3597482e-18/(hbar*omega)
    E_mol_upper[i] = E_mol[4]*4.3597482e-18/(hbar*omega)
    
    a[i] = BuschFunc.a_B_97_1st_branch(B[i], V_L)
    

    # Jose calc -> [:,0] - r in rvDw, [:,1] - r*psi
    wf = np.loadtxt("/Users/UniversityThings/Thywissen Lab/MSc Project/Data-40K-ab/who="+("%.2f" % freq)+"kHz/Bfield="+("%.7f" % B[i])+"/fort."+str(fort_num))
    wf_upper = np.loadtxt("/Users/UniversityThings/Thywissen Lab/MSc Project/Data-40K-ab/who="+("%.2f" % freq)+"kHz/Bfield="+("%.7f" % B[i])+"/fort."+str(fort_num_upper))
    
    r = wf[:,0]*65.0223142660582
    cc_norm[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[:,2])**2,x=r))
    oc_norm[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[:,1])**2,x=r))
    oc_norm_shortrange[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[0:2507,1])**2,x=r[0:2507]))
    oc_norm_shortrange[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[0:2519,1])**2,x=r[0:2519]))
    oc_norm_turningpoint[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[2518:2521,1])**2,x=r[2518:2521]))
    cc_norm_turningpoint[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf[2518:2521,2])**2,x=r[2518:2521]))
    
    cc_norm_upper[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf_upper[:,2])**2,x=r))
    oc_norm_upper[i] = np.sqrt(4*np.pi*scinteg.trapezoid(np.abs(wf_upper[:,1])**2,x=r))
    
    
CC_contact_lower = -8*np.pi*np.gradient(E_mol_lower,1/a)/a_ho
CC_contact_upper = -8*np.pi*np.gradient(E_mol_upper,1/a)/a_ho

CC_Nmol_lower = N_mol_func(B,CC_contact_lower)
CC_Nmol_upper = N_mol_func(B,CC_contact_upper)
    
fig = plt.figure(dpi=350,figsize=(8,6))

scale_CC = 3500
scale_CC_tp = 4e11
scale_OC_tp = 3500
scale_OC_sr = 3500

#plt.plot(B,cc_norm)
#plt.plot(B,oc_norm)
#plt.plot(B,C_upper)
#plt.plot(B,N_mol_upper)

fig= plt.figure(dpi=500,figsize=(6,6))
gs = gridspec.GridSpec(5,1) 

ax = plt.subplot(gs[0:3,0])
ax.plot(B,(cc_norm/np.sqrt(cc_norm**2+oc_norm**2))**2*scale_CC,linewidth=linewidth,label=r'CC $ar$ amp.')
#ax.plot(B,(oc_norm/np.sqrt(cc_norm**2+oc_norm**2))**2*scale_CC/100,label=r'CC $ab$ rel. amp.')
ax.plot(B,(oc_norm_shortrange/np.sqrt(cc_norm**2+oc_norm**2))**2*scale_OC_sr,'g-',linewidth=linewidth,label=r'CC $ab$ s.r. amp.')
ax.plot(B,(oc_norm_turningpoint/np.sqrt(cc_norm**2+oc_norm**2))**2*scale_OC_tp,'r:',linewidth=linewidth,label=r'CC $ab$ t.p. amp.')
ax.plot(B,(cc_norm_turningpoint/np.sqrt(cc_norm**2+oc_norm**2))**2*scale_CC_tp,'b:',linewidth=linewidth,label=r'CC $ar$ t.p. amp.')
ax.plot(B_lower[:-15000],N_mol[:-15000]*scale_CC,'--',linewidth=linewidth,label='PsP' )
#ax.plot(B_lower_no_eff_r[:-15000],N_mol_no_eff_r[:-15000]*scale_CC,'--',label='PsP w/ out eff. range')

#plt.plot(B_lower,N_mol,'.')

plot_data = 1

if plot_data ==1:
    
    #### 11 25 data ####
    P_max = 1
    
    power = 1.16
    data = scio.loadmat('/Volumes/GoogleDrive/My Drive/Lattice Shared/SharedData/2022 PA experiment/lattice_lifetime_all/Lattice_Lifetime_1mW_200_Er.mat')
    data = data['summary0']
    B_data = data[0][0][0][0,:] + 0.15 #G
    gamma_data = (data[0][0][1][0])*P_max/power #1/ms
    gamma_err = (data[0][0][2][0])*P_max/power #1/ms
    lifetimes = 1/gamma_data #ms
    lifetimes_err = gamma_err/gamma_data**2 #ms
    
    ax.errorbar(B_data,gamma_data,yerr=gamma_err,marker='o',linestyle='none',label='1.16 mW')
    
    #### 12 07 data ####
    power = 0.027
    data = scio.loadmat('/Volumes/GoogleDrive/My Drive/Lattice Shared/SharedData/2022 PA experiment/lattice_lifetime_all/Lattice_Lifetime_30uW_200_Er.mat')
    data = data['summary1']
    B_data07 = data[0][0][0][0,:]+0.15 #G
    gamma_data07 = (data[0][0][1][0])*P_max/power #1/ms
    gamma_err07 = (data[0][0][2][0])*P_max/power #1/ms
    
    ax.errorbar(B_data07,gamma_data07,yerr=gamma_err07,marker='o',linestyle='none',label='0.027 mW')


B_data_all = np.concatenate([B_data,B_data07])
gamma_data_all = np.concatenate([gamma_data,gamma_data07])
gamma_err_all = np.concatenate([gamma_err,gamma_err07])

#ax.set_xlim([195,208.5])
ax.set_xlim([202,208.5])
ax.set_ylim([0.5,500])


ax.set_xlabel('B field (G)')
ax.set_ylabel(r'Normalized photoexcitation rate (ms$^{-1}$mW$^{-1}$)')

#ax.legend(loc='lower left')
ax.legend(loc='upper right')


ax.set_yscale('log')

#residuals
data_ind = [68,74,51,57,62]
data_ind07 = [62,57,51,48,45,68,42]
data_ind_all = [68,74,51,57,62,62,57,51,48,45,68,42]
ax1 = plt.subplot(gs[3:5,:])

ax1.errorbar(B_data_all,gamma_data_all-(cc_norm[data_ind_all]/np.sqrt(cc_norm[data_ind_all]**2+oc_norm[data_ind_all]**2))**2*scale_CC,yerr=gamma_err_all,marker='o',linestyle='none',label=r'CC $ar$ amp.')
ax1.errorbar(B_data_all,gamma_data_all-(oc_norm_shortrange[data_ind_all]/np.sqrt(cc_norm[data_ind_all]**2+oc_norm[data_ind_all]**2))**2*scale_OC_sr,yerr=gamma_err_all,marker='o',linestyle='none',elinewidth=2,mec='g',mfc='g',ecolor='g',label=r'CC $ab$ s.r. amp.')
#ax1.errorbar(B_data_all,gamma_data_all-(oc_norm_turningpoint[data_ind_all]/np.sqrt(cc_norm[data_ind_all]**2+oc_norm[data_ind_all]**2))**2*scale_OC_tp,yerr=gamma_err_all,marker='o',linestyle='none',mec='r',mfc='r',ecolor='r',label=r'CC $ab$ t.p. amp.')
#ax1.errorbar(B_data_all,gamma_data_all-(cc_norm_turningpoint[data_ind_all]/np.sqrt(cc_norm[data_ind_all]**2+oc_norm[data_ind_all]**2))**2*scale_CC_tp,yerr=gamma_err_all,marker='o',linestyle='none',mec='b',mfc='b',ecolor='b',label=r'CC $ar$ t.p. amp.')

ax1.set_xlim([202,208.5])
ax1.set_ylim([-10,10])
#ax1.set_yscale('log')

ax1.set_xlabel('B field (G)')
ax1.set_ylabel(r'Residuals (ms$^{-1}$mW$^{-1}$)')

plt.tight_layout()

#fig.suptitle(r'$\omega_{\mathrm{HO}} = 2\pi \times 128.66$ kHz lower branch ')

#Contact figure
fig, ax = plt.subplots(dpi=500,figsize=(6,6))
ax.plot(202.15*np.ones(10),np.linspace(0,100,10),'--',linewidth=1,color='limegreen',label='Feshbach resonance')

ax.plot(B,CC_contact_lower*a_ho,'.',mfc='mediumslateblue',mec='mediumslateblue',label='B0: CC adiab.')
ax.plot(B_lower[0:58100],C_array[0:58100,2],'--',color='blueviolet',linewidth=2,label='B0: PsP')

ax.plot(B,CC_contact_upper*a_ho,'.',mfc='darkorange',mec='darkorange',label='B1: CC adiab.')
ax.plot(B_upper[2000:87000],C_array_upper[2000:87000,2],'--',color='maroon',linewidth=2,label='B1: PsP')

#ax.plot(B,BuschFunc.contact_lower_interp(B, V_L))

#ax.set_xlim([195,210])
ax.set_xlim([196,208])
ax.set_ylim([0,65])
#ax.set_ylim([0,0.020])
ax.set_xlabel('B field (G)')
ax.set_ylabel(r'Contact ($a_{\mathrm{ho}}^{-1}$)')

ax.legend(loc='upper right')

plt.tight_layout(pad=2)

#N mol figure
N_mol_interp_func = scinterp.interp1d(B_lower,N_mol,'linear')
N_mol_interp = N_mol_interp_func(B)

N_mol_no_eff_r_interp_func = scinterp.interp1d(B_lower_no_eff_r,N_mol_no_eff_r,'linear')
N_mol_no_eff_r_interp = N_mol_no_eff_r_interp_func(B)

N_mol_upper_interp_func = scinterp.interp1d(B_upper,N_mol_upper,'linear')
N_mol_upper_interp = N_mol_upper_interp_func(B)

N_mol_no_eff_r_upper_interp_func = scinterp.interp1d(B_upper_no_eff_r,N_mol_no_eff_r_upper,'linear')
N_mol_no_eff_r_upper_interp = N_mol_no_eff_r_upper_interp_func(B)


fig = plt.figure(dpi=500,figsize=(6,6))
gs = gridspec.GridSpec(3,1) 

ax = plt.subplot(gs[0:2,0])


ax.plot(202.15*np.ones(10),np.linspace(0,100,10),'--',linewidth=1,color='limegreen',label='Feshbach resonance')

ax.plot(B,(cc_norm/np.sqrt(cc_norm**2+oc_norm**2))**2,linewidth=4,color='mediumblue',label=r'B0: CC w.f. amp.')
ax.plot(B,CC_Nmol_lower,'.',mfc='mediumslateblue',mec='mediumslateblue',label='B0: CC adiab.')



ax.plot(B_lower[:-15000],N_mol[:-15000],'--',color='blueviolet',linewidth=2,label='B0: PsP' )
ax.plot(B_lower_no_eff_r[:-15000],N_mol_no_eff_r[:-15000],'--',linewidth=2,color='deepskyblue',label='B0: PsP no eff. range' )

ax.plot(B,(cc_norm_upper/np.sqrt(cc_norm_upper**2+oc_norm_upper**2))**2,linewidth=4,color='red',label=r'B1: CC w.f. amp.')
ax.plot(B,CC_Nmol_upper,'.',mfc='darkorange',mec='darkorange',label='B1: CC adiab.')
ax.plot(B_upper[:-15000],N_mol_upper[:-15000],'--',color='maroon',linewidth=2,label='B1: PsP' )

ax.plot(B_upper_no_eff_r[:-15000],N_mol_no_eff_r_upper[:-15000],'--',linewidth=2,color='yellow',label='B1: PsP no eff. range' )

ax.set_xlim([196,208])
ax.set_ylim([0,0.12])

ax.set_xlabel('B field (G)')
ax.set_ylabel(r'$Z$ (arb.)')

ax.legend(loc='upper right')



ax1 = plt.subplot(gs[2,0])

ax1.plot(B,(CC_Nmol_lower-(cc_norm/np.sqrt(cc_norm**2+oc_norm**2))**2)/CC_Nmol_lower*100,linewidth=4,color='mediumblue',label=r'B0: CC w.f. amp.')
ax1.plot(B,(CC_Nmol_lower-N_mol_interp)/CC_Nmol_lower*100,'--',color='blueviolet',linewidth=2,label='B0: PsP' )
ax1.plot(B,(CC_Nmol_lower-N_mol_no_eff_r_interp)/CC_Nmol_lower*100,'--',linewidth=2,color='deepskyblue',label='B0: PsP no eff. range' )

ax1.plot(B,(CC_Nmol_upper-(cc_norm_upper/np.sqrt(cc_norm_upper**2+oc_norm_upper**2))**2)/CC_Nmol_upper*100,linewidth=2,color='red',label=r'B1: CC w.f. amp.')
ax1.plot(B,(CC_Nmol_upper-N_mol_upper_interp)/CC_Nmol_upper*100,'--',color='maroon',linewidth=2,label='B1: PsP' )
ax1.plot(B,(CC_Nmol_upper-N_mol_no_eff_r_upper_interp)/CC_Nmol_upper*100,'--',linewidth=2,color='yellow',label='B1: PsP no eff. range' )

ax1.set_xlabel('B field (G)')
ax1.set_ylabel(r'% diff. from $Z_{\mathrm{CC \ adiab.}}$')


ax1.set_xlim([196,208])
ax1.set_ylim([-10,10])
plt.tight_layout(pad=2)



    
    
