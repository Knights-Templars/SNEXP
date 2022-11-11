#!/usr/bin/env python
# coding: utf-8

# In[95]:


import os


# In[96]:


version = 1.0
model_file = "sed_model_" + str(version) + "_1.7e+51.mod"
param_file = "param_model_" + str(version) + "_1.7e+51.lua"


# In[97]:


folder_path = '/Users/anirbandutta/Documents/SNEXP/Violent_Merger_2012_Model_Ejecta/'


# In[98]:


f = open(folder_path + param_file, 'w')
f.write('-- param file for running the light curve of a 1D Type Ia-like supernova' + '\n')
f.write('-- using LTE and line expansion opacity' + '\n')
f.write('-- atomic line data taken from kurucz list ("fuzz")' + '\n')
f.write("\n")
f.write('-- model type and file' + '\n')
f.write('grid_type = "grid_1D_sphere"' + '\n')
f.write('model_file = "%s"'%model_file + '\n')
f.write('hydro_module = "homologous"' + '\n')
f.write("\n")
f.write('-- home directory and atomic data files' + '\n')
f.write("\n")
f.write('sedona_home = os.getenv("SEDONA_HOME")' + '\n')
f.write('defaults_file = sedona_home.."/defaults/sedona_defaults.lua"' + '\n')
f.write('data_atomic_file   = sedona_home.."/data/ASD_atomdata.hdf5"' + '\n')
f.write('data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"' + '\n')
f.write("\n")
f.write("-- helper variable" + '\n')
f.write("days = 3600.0*24" + '\n')
f.write("\n")
f.write("-- total number of particles used initially" + '\n')
f.write("particles_n_initialize = 1e6" + '\n')
f.write("-- number of particles emitted per time step from radioactivity" + '\n')
f.write("particles_n_emit_radioactive = 1e6" + '\n')
f.write("particles_last_iter_pump = 10" + '\n')
f.write("\n")
f.write("-- time start/stop and stepping" + '\n')
f.write("tstep_time_start = 1*days" + '\n')
f.write("tstep_time_stop = 100*days" + '\n')
f.write("tstep_max_dt = 1.0*days" + '\n')
f.write("tstep_max_delta = 0.1" + '\n')
f.write("\n")
f.write("-- frequency grid to calculate and store opacities" + '\n')
f.write("nu1 = 1e14" + '\n')
f.write("nu2 = 2e16" + '\n')
f.write("transport_nu_grid = {nu1,nu2,0.0003,1}" + '\n')
f.write("\n")
f.write("-- frequency grid to calculate output spectrum" + '\n')
f.write("nu1_spec = nu1*1.1" + '\n')
f.write("spectrum_nu_grid = {nu1_spec,nu2,0.002,1}" + '\n')
f.write("spectrum_time_grid = {0*days,100*days,1.0*days}" + '\n')
f.write("output_time_write = 1*days" + '\n')
f.write("output_write_radiation = 1" + '\n')
f.write("\n")
f.write("-- opacity settings" + '\n')
f.write("opacity_grey_opacity = 0.0" + '\n')
f.write("opacity_epsilon = 1.0" + '\n')
f.write("opacity_electron_scattering = 1" + '\n')
f.write("opacity_free_free = 1" + '\n')
f.write("opacity_bound_bound = 1" + '\n')
f.write("opacity_bound_free = 1" + '\n')
f.write("opacity_line_expansion = 1" + '\n')
f.write("opacity_fuzz_expansion = 1" + '\n')
f.write("\n")
f.write("-- transport settings" + '\n')
f.write("transport_steady_iterate = 0" + '\n')
f.write("transport_radiative_equilibrium = 1" + '\n')


f.close()


# In[ ]:





# In[ ]:




