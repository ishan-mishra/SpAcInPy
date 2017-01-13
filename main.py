#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 13:24:02 2016

@author: ishan
"""

'''StAcInPy'''

'''Tool for calculating spectral indices for a given stellar spectra'''

import configparser
from astropy.io import fits, ascii
from astropy.constants import c
from lib import *
import numpy as np
import sys

#import necessary information from the configuration file param.ini
#param.ini should be in the same directory as the SpAcInPy script

config = configparser.ConfigParser()

config.read('param.ini')

#Find out the input file type (e2ds, s1d or ascii table)

file_type = str(config['FILE TYPE']['file_type'])

print '\n Input file type:', file_type, '\n'

if file_type == 'e2ds':
    
    wav_path = str(config['E2DS INPUTS']['wav_file'])
    ccf_path = str(config['E2DS INPUTS']['ccf_file'])
    e2ds_path = str(config['E2DS INPUTS']['e2ds_file'])
    

    #Information about the RVs from the header of the ccf fits file

    star_rv = str(config['E2DS INPUTS']['star_rv'])
    earth_rv = str(config['E2DS INPUTS']['earth_rv'])
    
    #Check whether the keys are floats (actual RVs) or strings (keywords for header)
    
    try:
        RVC = float(star_rv)
        BERV = float(earth_rv)
    
    except:
        pass
    
    try:
        hdu_ccf = fits.open(ccf_path)
        RVC = hdu_ccf[0].header[star_rv]
        BERV = hdu_ccf[0].header[earth_rv]

        hdu_ccf.close()
    
    except ValueError:
        
        print 'Unacceptable inputs for RV values'
        

    #The total RV (in km/s) is the sum of the two RVs obtained

    RV_tot = RVC + BERV 

    #Read in the wav file and correct all the pixel values (wavelengths) using RV_tot and
    #the formula for Doppler shift.

    hdu_wav = fits.open(wav_path)

    wav_soln = hdu_wav[0].data
 
    wav_soln_corr = wav_soln + (wav_soln*RV_tot*1000/c.value)

    hdu_wav.close()


    #The spectra

    hdu_e2ds = fits.open(e2ds_path)

    spectra = hdu_e2ds[0].data

    #Information about the spectral lines

    line_list = [str(item) for item in config['LINES']]
    line_dicts = getLineDict(line_list,config,'LINES')

    #Locate the wavelength range in terms of pixels (region of interest around the 
    #center wavelength) in the corrected  wavelength solution. Inform whether the 
    #wavelength range is present in just 1 order or 2 consecutive orders.


    for line in line_dicts:
    
        lower = line['center'] - line['roi']/2.0
        upper = line['center'] + line['roi']/2.0

        #ensure that given wavelenghts are within the bounds of the soln
    
        if lower < wav_soln_corr.min() or upper > wav_soln_corr.max():
           sys.exit('Wavelength limits is out of bounds!')
    
        #Add a new entry in the line's dictionary for the obtained pixel range 
        line['pixel_range'] = get2DPixRange(line['name'],lower,upper,wav_soln_corr)
     

    for item in line_dicts:
        print item, '\n'
    
        
    #Calculating the spectral index. We are calculating the R_HK index here, using
    #the formulae and calibrations described in the 2011 paper titled 'HARPS search for 
    #southern extra-solar planets' by Lovis et al.
    
    
    #Get the passband locations

    #Information about the passband regions

    pb_list = [str(item) for item in config['PASSBANDS']]
    pb_dicts = getLineDict(pb_list,config,'PASSBANDS')
    
    #Locate the wavelength range in terms of pixels.
    
    for pb in pb_dicts:
    
        lower = pb['center'] - pb['roi']/2.0
        upper = pb['center'] + pb['roi']/2.0

        #ensure that given wavelenghts are within the bounds of the soln
    
        if lower < wav_soln_corr.min() or upper > wav_soln_corr.max():
           sys.exit('Wavelength limits is out of bounds!')
    
        #Add a new entry in the line's dictionary for the obtained pixel range 
        pb['pixel_range'] = get2DPixRange(pb['name'],lower,upper,wav_soln_corr)    
        
    
    #Calculate the log(R_HK) value
    
    #Calculate flux in H,K,R and V
    
    K_flux = sumFlux(line_dicts[0],wav_soln_corr, spectra)
    H_flux = sumFlux(line_dicts[1],wav_soln_corr, spectra)
    V_flux = sumFlux(pb_dicts[0],wav_soln_corr, spectra)
    R_flux = sumFlux(pb_dicts[1], wav_soln_corr, spectra)
    
    
    
    b_v = float(str(config['MISC']['b_v']))
    
    #If two different total flux values are present for a passband, chose the larger one
    
    
    R_HK = calcR_HK(max(H_flux),max(K_flux),max(V_flux),max(R_flux), b_v)
    logR_HK = np.log10(R_HK)
    
    print R_HK, logR_HK 
    
    
    #Calculate H_alpha index value
    
    H_alpha_flux = sumFlux(line_dicts[2],wav_soln_corr, spectra)
    L_alpha_flux = sumFlux(pb_dicts[2],wav_soln_corr, spectra)
    R_alpha_flux = sumFlux(pb_dicts[3],wav_soln_corr, spectra)
    
    H_alpha_index = max(H_alpha_flux)/(max(L_alpha_flux) + max(R_alpha_flux))
            
    
    hdu_e2ds.close()
    
 
elif file_type == 's1d':
    
    s1d_path = str(config['S1D INPUTS']['s1d_file'])
    hdu_s1d = fits.open(s1d_path)
    
    #Information about the RVs from the header of the ccf fits file

    #star_rv = str(config['S1D INPUTS']['star_rv'])
    earth_rv = str(config['S1D INPUTS']['earth_rv'])
    
    #Check whether the keys are floats (actual RVs) or strings (keywords for header)
    
    try:
        #RVC = float(star_rv)
        BERV = float(earth_rv)
    
    except:
        pass
    
    try:
        #RVC = hdu_s1d[0].header[star_rv]
        BERV = hdu_s1d[0].header[earth_rv]
    
    except ValueError:
        
        print 'Unacceptable inputs for RV values'
        

    #The total RV (in km/s) is the sum of the two RVs obtained

    #RV_tot = RVC + BERV 
    RV_tot = BERV
    
    
    #Information about the spectral lines

    line_list = [str(item) for item in config['LINES']]
    line_dicts = getLineDict(line_list,config)
        
    #The spectra

    spectra = hdu_s1d[0].data
    

    #Define the wavelength axis
    
    #Extract wavelength parameters
    
    wav_ref = (str(config['S1D INPUTS']['wav_ref']))
    wav_step = (str(config['S1D INPUTS']['wav_step']))
    
    #Check whether the keys are floats (actual wavelength parameters) or strings
    #(keywords for header)
    
    try:
        wav_ref = float(wav_ref)
        wav_step = float(wav_step)
    
    except:
        pass
    
    try:
        wav_ref = hdu_s1d[0].header[wav_ref]
        wav_step = hdu_s1d[0].header[wav_step]
    
    except ValueError:
        
        print 'Unacceptable inputs for wavelength parameters'
        
    #Define wavelength axis
      
    wav_axis = np.empty_like(spectra)

    wav_axis[0] = wav_ref

    for i in range(1,len(wav_axis)):
        
        wav_axis[i] = wav_axis[i-1] + wav_step

    #Find locations of the spectral lines
    
    for line in line_dicts:
    
        lower = line['center'] - line['roi']/2.0
        upper = line['center'] + line['roi']/2.0

        #ensure that given wavelenghts are within the bounds of the wav_axis
    
        if lower < wav_axis.min() or upper > wav_axis.max():
            sys.exit('Wavelength limit is out of bounds!')
    
        #Add a new entry in the line's dictionary for the obtained pixel range 
        line['pixel_range'] = get1DPixRange(line['name'],lower,upper,wav_axis)
     

    for item in line_dicts:
        print item, '\n'    
        
        
    #Calculating the spectral index. We are calculating the R_HK index here, using
    #the formulae and calibrations described in the 2011 paper titled 'HARPS search for 
    #southern extra-solar planets' by Lovis et al.
    
    
    #Get the passband locations
        
    
    hdu_s1d.close()

    
elif file_type == 'ascii':
    
    ascii_path = str(config['ASCII INPUTS']['ascii_file'])
    data = ascii.read(ascii_path)
    cols = data.colnames
    
    wav_axis = data[cols[0]]
    
    #Convert wavelength units to angstrom
    
    wav_axis = wav_axis*10000.0
    flux = data[cols[1]]

    #Information about the spectral lines

    line_list = [str(item) for item in config['LINES']]
    line_dicts = getLineDict(line_list,config)    
    
    #Find locations of the spectral lines
    
    
    for line in line_dicts:
    
        lower = line['center'] - line['roi']/2.0
        upper = line['center'] + line['roi']/2.0
        
        #ensure that given wavelenghts are within the bounds of the wav_axis
    
        if lower < wav_axis.min() or upper > wav_axis.max():
            sys.exit('Wavelength limits is out of bounds!')    
    
        #Add a new entry in the line's dictionary for the obtained pixel range 
        line['pixel_range'] = get1DPixRange(line['name'],lower,upper,wav_axis)
     

    for item in line_dicts:
        print item, '\n'    