#!/usr/bin/env python
# coding: utf-8

# In[79]:


import numpy as np
import pandas as pd
from astropy.constants import M_sun
import matplotlib.pyplot as plt

plot = True
save = False
mass = 1.0        # mass in Msolar
KE = 1.7e51       # KE in erg

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['axes.titlesize'] = 25
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['legend.fontsize'] = 20
plt.rc('font', family='Times New Roman')

nx = 100          # number of radial zones
texp = 86400      # secs

ve = (KE / (6 * mass * M_sun.value * 10**3))**0.5                                      # cm/s
rho_0 = (mass * M_sun.value * 10**3) / (8 * np.pi * (ve * texp)**3.0)       # gm/cm^3

v_max = 4.0e9     # outer velocity (cm/s)
dv = v_max / (1.0 * nx)

rho = np.zeros(nx)
v = np.zeros(nx)

for i in range(nx):
    
    v[i] = (i + 1.0) * dv
    vm = (i + 0.5) * dv
    
    rho[i] = rho_0 * np.exp(-vm / ve)

model_df = pd.DataFrame({'velocity': v, 'density': rho})

if plot:
	model_df.plot('velocity', 'density')
	plt.yscale("log")
	plt.xscale("log")

plt.show();

if save:
    model_df.to_csv('/Users/anirbandutta/Documents/SNEXP/Violent_Merger_2012_Model_Ejecta/' + 'model_' + str(mass) + '_' + 
                 str(KE), sep='\t')





