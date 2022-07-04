#!/usr/bin/env python
# coding: utf-8

# This script will calculate airmass and put it in header.
# This is a standalone script to calculate airmass besides spec_reduce_v1.py

#Author: Anirban Dutta
#Version: v1.0

#===========================================================================================#

import os
import re
import glob
import datetime
import warnings
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
warnings.filterwarnings("ignore")
#===========================================================================================#

# Change current working directory and object details 
# This should be provided by the user
# the filename should have OBJECT NAME to run this code, e.g. 2020-10-15_SN2020uxz.fits
#===========================================================================================#
cwd = ''
OBJECT_NAME='feige34'           # Name of the object
OBJECT_RA = '10:39:36.71'       # Object RA
OBJECT_DEC = '+43:06:10.1'      # Object DEC
#===========================================================================================#


# HFOSC old CCD specifications

read_noise=4.87
gain=1.22
data_max=55000

# HFOSC new CCD specifications

read_noise_new = 5.75
gain_new = 0.28
data_max_new = 700000

#image header keywords
# Add keywords here which you want to see

Date = 'DATE-AVG'
Filter = 'FILTER'
OBJECT = 'OBJECT'
CCDSECTION = 'CCDSEC'
RA = 'RA'
DEC = 'DEC'



# Observatory Details
#===========================================================================================#


observatory = 'iao'
obs_lat = '32:46:46'
obs_long = '78:57:51'
latitude = 32.7794 * u.deg
longitude = 78.9642 * u.deg
altitude = 4500 * u.m
tz = +5.5 * u.hour
name = 'Indian Astronomical Observatory, Hanle'

#===========================================================================================#


# In[6]:


Software_Code = 'spec_reduce1_v1.0'
Author_Name = 'Anirban Dutta'

#===========================================================================================#


def remove_file(file_name):
	
	try:
		os.remove(file_name)
	except OSError:
	    pass

# Function for removing files having similar names: argument: common_text

def remove_similar_files(common_text):
	
	for residual_file in glob.glob(common_text):
		remove_file(residual_file)
		
		

# Function to group similar files: arguments: text_list, common_text, exceptions

def group_similar_files(text_list, common_text, exceptions = ''):
    
    '''
    text_list: A text file used to store list of files
    common_text: A string (e.g. *.fits, *.list) used for grouping similar
    kinds of files
    exceptions: string of file name to exclude in grouping
    
    returns: list of grouped files
    '''
    
    list_files = glob.glob(common_text)
    if exceptions != '':
    	list_exceptions = exceptions.split(',')
    	for text in list_exceptions:
        	list_files = filter(lambda x: not re.search(text, x), list_files)
        
    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files        


def text_list_to_python_list(text_list):
	if os.path.exists(text_list):
		with open(text_list, 'r+') as f:
			python_list=f.read().split()

	return python_list
	
def python_list_to_text_list(python_list, text_list):
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')
            
def list_lists_to_list(list_lists, text_list):

    list_name=[]
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list=f.read().split()
            for element in file_list:
                list_name.append(element)
    python_list_to_text_list(list_name, text_list)
    
    return list_name
    
#===========================================================================================#

def copy_header(file_name):
        
    headerlist = fits.open(file_name, mode = 'update')        
    file_header = headerlist[0].header
    date_obs = file_header['DATE-OBS']
    exposure_time = file_header['EXPTIME']
    Filter = file_header['FILTER']
    
    object_ra = RA
    object_dec = DEC
    Object = OBJECT_NAME
            
        
    if re.search(':', object_ra):
        c = SkyCoord(object_ra, object_dec, unit=(u.hourangle, u.deg))
        ra_deg = c.ra.degree
        dec_deg = c.dec.degree
        RA_degrees = round(ra_deg, 5)
        DEC_degrees = round(dec_deg, 5)
    else:
        print("The RA and DEC are already in degrees") 
        
        RA_degrees == object_ra 
        DEC_degrees == object_dec     
        
                
    list_keywords = ['DATE-OBS', 'EXPTIME', 'FILTER', 'RA', 'DEC', 'OBJECT', 'RA_DEG', 'DEC_DEG'] 
    dict_header = {'DATE-OBS': date_obs, 'EXPTIME': exposure_time, 'FILTER': Filter, 'RA': object_ra, 'DEC': object_dec,                  'OBJECT': Object, 'RA_DEG': RA_degrees, 'DEC_DEG': DEC_degrees}
    comment = {'DATE-OBS': 'The date of observation', 'EXPTIME': 'Time of exposure of the frame', 'FILTER': 'Filter',                 'RA': 'Right Ascension of the target', 'DEC': 'Declination of the target', 'OBJECT': 'Name of the Object under study',                 'RA_DEG': 'Right Ascension in degrees', 'DEC_DEG': 'Declination in degrees'}
    
    for keyword in list_keywords:
        if keyword in file_header.keys():
            file_header.remove(keyword, remove_all = True)
        file_header.append(card = (keyword, dict_header[keyword], comment[keyword])) 
        
    headerlist.flush()
    headerlist.close()                 

def edit_header(textlist_files, object_name):
    
    list_files = text_list_to_python_list(textlist_files)
    
    object_ra = RA
    object_dec = DEC
    Object = OBJECT_NAME
    
    for file_name in list_files:
        hdulist = fits.open(file_name, mode = 'update')
        file_header = hdulist[0].header
        OBJECT = file_header['OBJECT']
        
        
        
        list_keywords = ['OBJECT', 'RA', 'DEC']
        dict_header = {'OBJECT': OBJECT, 'RA': RA, 'DEC': DEC}
       
        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all = True)

        list_update_keywords = ['OBJECT', 'RA', 'DEC', 'COMMENT', 'AUTHOR']
        dict_update_header = {'OBJECT': object_name, 'RA': OBJECT_RA, 
                              'DEC': OBJECT_DEC, 'COMMENT': Software_Code, 'AUTHOR': Author_Name}

        for keyword in list_update_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all = True)
            file_header.append(card=(keyword, dict_update_header[keyword]))
        hdulist.flush()
        hdulist.close() 
        
        
def calculate_airmass(textlist_files):
    
    '''
    textlist_files: List of fits files in a text for which you want airmass to be put in the header
    
    returns: list of airmasses.
    Works for both old and new ccd on HFOSC.
    '''
    
    list_airmass = []
    
    list_files = text_list_to_python_list(textlist_files)
    hct = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude)
    for file_name in list_files:
        hdu = fits.open(file_name, mode='update')
        header = hdu[0].header

        if 'TM_START' in header.keys():
            date_obs = header['DATE-OBS']
            time_start = header['TM_START']
            
            ra_ = OBJECT_RA
            dec_ = OBJECT_DEC
            c = SkyCoord(ra_, dec_, unit=(u.hourangle, u.deg))
            ra = c.ra.deg
            dec = c.dec.deg
            time_utc = str(datetime.timedelta(seconds=int(time_start)))
            datetime_utc = date_obs+' '+time_utc
            time = Time(datetime_utc) 
        else:
            date_time = header['DATE-AVG'].split('T')
            time_obj = date_time[0]+' '+date_time[1] 
            time = Time(time_obj)
            
            ra = OBJECT_RA
            dec = OBJECT_DEC
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            ra = c.ra.deg
            dec = c.dec.deg
        

        coord = SkyCoord(ra, dec, unit='deg')
        altaz_ = coord.transform_to(AltAz(obstime=time, location=hct))
        airmass = altaz_.secz.value
        print('The image %s has been observed at an airmass of %f'%(file_name, airmass))
        
        list_keywords = ['COMMENT', 'RA', 'DEC', 'AIRMASS']
        dict_header = {'COMMENT':Software_Code, 'RA':ra, 'DEC':dec,
                      'AIRMASS': airmass}
        
        for key in list_keywords:
            if key in header.keys():
                header.remove(key, remove_all=True)
            header.append(card=(key, dict_header[key]))
            
        hdu.flush()
        hdu.close()
        
        list_airmass.append(airmass)
        
    return list_airmass   

#===========================================================================================#

print('The current working directory is:\n %s' % cwd)

os.chdir(cwd)

#===========================================================================================#

list_object = group_similar_files('list_object_spec',  '*'+OBJECT_NAME+'*.fits')
edit_header('list_object_spec', OBJECT_NAME)
list_spec = group_similar_files('list_clean_spec', '*'+OBJECT_NAME+'*.fits')
list_airmass = calculate_airmass(textlist_files='list_clean_spec')

#===========================================================================================#

