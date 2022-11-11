#!/usr/bin/env python
# coding: utf-8

# In[85]:


import os
import re
import glob
import shutil

working_directory = os.getcwd()
os.chdir(working_directory)
print(working_directory)

list_param_files = glob.glob(working_directory + '/' + '*.lua')
sedona_executable = 'sedona6.ex'

list_model_files = glob.glob(working_directory + '/' + '*.mod')


for i in range(len(list_model_files)):
    
    file = list_model_files[i].split('/')[-1]
    name = file[4:21]
    
    for j in range(len(list_param_files)):
        
        param_file = list_param_files[j].split('/')[-1]
        if name in param_file:
            if os.path.exists(name):
                shutil.rmtree(name)
            os.mkdir(name)
            shutil.copy(file, name)
            shutil.copy(param_file, name)
            shutil.copy(sedona_executable, name)
        
 
for root, subdirectories, files in os.walk(working_directory):
    for i in subdirectories:
        os.chdir(root + '/' + i)
        for j in glob.glob('param_*'):
            print ("Running Sedona in %s"%os.getcwd())
            command = "mpirun -n 8 ./sedona6.ex %s"%j
            os.system(command)
            #print (os.listdir())

