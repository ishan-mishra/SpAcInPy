#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 11:41:54 2016

@author: ishan
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


def findNearest(array,value):
    
    '''Function finds the pixel containing the value closest to the input
       and returns its index'''
       
    idx = (np.abs(array-value)).argmin()
    if array.ndim == 1:
        return idx
    else: return (idx//array.shape[1],idx%array.shape[1])
    
    
    

def plotSpectra(spectra_file,wav_soln,label,start,end):
    
    '''Function to plot the pixel range for a given Spectra'''
    
    ''' Note: start and end are pixel indices'''
    
    idx_start = start
    idx_end   = end
    
    wav_array = wav_soln.flatten()
    spec_array = spectra_file.flatten()
    
    i_idx = idx_start[0]*wav_soln.shape[1] + idx_start[1]
    f_idx = idx_end[0]*wav_soln.shape[1] + idx_end[1]

    x = wav_array[i_idx:f_idx]
    y = spec_array[i_idx:f_idx]

    plt.plot(x,y)
    plt.title(label)
    
    
def returnFlat(array,start,end):
    
    '''To return a 1D array corresponding to elements of the input array between
       start and end pixels'''
    '''Note: start and end inputs are tuples'''
       
    flat = array.flatten()
    
    i_idx = start[0]*array.shape[1] + start[1]
    f_idx = end[0]*array.shape[1] + end[1]
    
    return flat[i_idx:f_idx]

def getLineDict(line_list,config,section):
    
    '''To retrun a list of dictionaries containing info about each spectral line
       correspoinding to the input list of spectral lines'''
       
    line_dicts = []
    
    for item in line_list:
        info = str(config[section][item]).split()
        line_dicts.append(dict([('name',info[0]),('center',float(info[1])), 
                             ('roi',float(info[2]))]))

    return line_dicts
    
    
def get2DPixRange(line_name,lower,upper,wav_soln_corr):
    
    '''To return the pixel locations corresponding to the 'lower' and 'upper' wavelength
       limits for a spectral line in a e2ds file. The exact wavelengths may not be present
       in the solution, so the pixels with the closest values are returned.'''
       
       
    #find the indices of the two pixel values closest to the input wavelenghts
    
    idx_l = [np.argsort(np.abs(wav_soln_corr.flatten() - lower))[0],
           np.argsort(np.abs(wav_soln_corr.flatten() - lower))[1]]
     
    idx_l[0] = (idx_l[0]//wav_soln_corr.shape[1],idx_l[0]%wav_soln_corr.shape[1])
    idx_l[1] = (idx_l[1]//wav_soln_corr.shape[1],idx_l[1]%wav_soln_corr.shape[1])
    
    idx_l = sorted(idx_l, key=lambda x: x[0])
    #same for the upper wavelength
    
    idx_u = [np.argsort(np.abs(wav_soln_corr.flatten() - upper))[0],
           np.argsort(np.abs(wav_soln_corr.flatten() - upper))[1]]
     
    idx_u[0] = (idx_u[0]//wav_soln_corr.shape[1],idx_u[0]%wav_soln_corr.shape[1])
    idx_u[1] = (idx_u[1]//wav_soln_corr.shape[1],idx_u[1]%wav_soln_corr.shape[1])    
    
    idx_u = sorted(idx_u, key=lambda x: x[0])
    
    #Check whether the wavelength range exists in two consecutive rows/orders
    
    if idx_l[0][0] != idx_l[1][0] and idx_u[0][0] != idx_u[1][0] and \
    idx_l[0][0] == idx_u[0][0] and idx_l[1][0] == idx_u[1][0]:
        
        print '\n Wavelength range for ',line_name,' exists in two different orders: \n'
        
        
        print 'First wavelength pixel range: ',(idx_l[0],idx_u[0])
        print 'Second wavelength pixel range: ',(idx_l[1],idx_u[1])
        
        return [(idx_l[0],idx_u[0]),(idx_l[1],idx_u[1])]
               
    else:
        
        print '\n Wavelength pixel range for ',line_name,'is present in one order: \n'
        
        #Ensure that lower wavelength and upper wavelength are not in different rows/orders
        if idx_l[0][0] == idx_u[0][0]:
            
            print (idx_l[0],idx_u[0]), '\n'
            
            return [(idx_l[0],idx_u[0])]
            
        elif idx_l[1][0] == idx_u[0][0]:
            
            print (idx_l[1],idx_u[0])
            return [(idx_l[1],idx_u[0])]
            
        else: print (idx_l[1],idx_u[1]), '\n'
        return [(idx_l[1],idx_u[1])]
            

def get1DPixRange(line_name,lower,upper,wav_axis):
    
    '''To return the pixel locations corresponding to the 'lower' and 'upper' wavelength
       limits for a spectral line in a e2ds file. The exact wavelengths may not be present
       in the solution, so the pixels with the closest values are returned.'''
     
       
    #find the indices of the two pixel values closest to the input wavelenghts
    
    idx_lower = (np.abs(wav_axis-lower)).argmin()
    idx_upper = (np.abs(wav_axis-upper)).argmin()

    return [(idx_lower, idx_upper)]    
                
                

def sumFlux(line_info,wav_soln_corr,spectra):
    
    '''Function to sum up the flux between a given wavelength range. 
       
       Info. about the spectral line center, the width of the region of interest around 
       the center,and the 'start' and 'end' pixel indices should be contained in the 
       'line_info' dictionary.'''
       
       
    #Create a list to store the calculated sum. A list is needed because the wavelength
    #range can be present in 2 orders.
    
    total_flux = []

    #Iterate over the 'pixel_range' entry in the dictionary and calculate flux for each
    #pixel range
    
    for pix_range in line_info['pixel_range']:
        
        start_pix = pix_range[0]
        end_pix = pix_range[1]

        #Calculate the exact wavelength limits using the spectral line info entered by user
        
        start_wav = line_info['center'] - line_info['roi']/2.0
        end_wav = line_info['center'] + line_info['roi']/2.0

        #Find the location of the start_wav and end_wav in the corrected wavelength
        #solution. Calculate the fraction of the pixels at the two ends of the range
        #that need to be included.
        
        #For start_wav
        
        if start_wav > wav_soln_corr[start_pix]:
            
            lambda_1 = wav_soln_corr[start_pix]
            lambda_2 = wav_soln_corr[start_pix[0],start_pix[1]+1] #the next wavelength
            
            f_start = (lambda_2 - start_wav)/(lambda_2 - lambda_1)
            
            #print 'f_start = ', f_start
            #print lambda_1, lambda_2, start_wav
            
        elif start_wav < wav_soln_corr[start_pix]:
            
            lambda_1 = wav_soln_corr[start_pix[0],start_pix[1]-1]
            lambda_2 = wav_soln_corr[start_pix]

            f_start = (lambda_2 - start_wav)/(lambda_2 - lambda_1)
            
            #print 'f_start = ', f_start
            
            #update the start_pix to one pixel before the current start_pix to include the 
            #fractional flux from this new start_pix
            start_pix = (start_pix[0],start_pix[1]-1)
            
        elif start_wav == wav_soln_corr[start_pix]:
            
            f_start = 1.0
            
            #print 'f_start = ', f_start
            
        #For end_wav
        
        if end_wav < wav_soln_corr[end_pix]:
            
            lambda_1 = wav_soln_corr[end_pix[0],end_pix[1]-1]
            lambda_2 = wav_soln_corr[end_pix] #the next wavelength
            
            f_end = (lambda_2 - end_wav)/(lambda_2 - lambda_1)
            
            #update the end_pix to one pixel before the current end_pix. The current end_pix
            #lies outside the wavelength range.
            
            end_pix = (end_pix[0],end_pix[1]-1) 
            
            #print 'f_end = ', f_end
            
        elif end_wav > wav_soln_corr[end_pix]:
            
            lambda_1 = wav_soln_corr[end_pix]
            lambda_2 = wav_soln_corr[end_pix[0],end_pix[0]]

            f_end = (lambda_2 - end_wav)/(lambda_2 - lambda_1)
            
            #print 'f_end = ', f_end
            
        elif end_wav == wav_soln_corr[end_pix]:
            
            f_end = 1.0        
        
            #print 'f_end = ', f_end
        
        
        #Now calculate the total flux in the wavelength range
        
        #Contribution of the intermediate pixels, ie, all the pixels except start_pix
        #and end_pix
        
        int_flux = 0.0
        
        for value in spectra[start_pix[0],start_pix[1]+1:end_pix[1]]:
            
            int_flux = int_flux + value
            
        #Contribution of start_pix
        
        start_flux = f_start*spectra[start_pix]

        #Contribution of end_pix
        
        end_flux = f_end*spectra[end_pix]


        #Total flux
        
        flux = start_flux + int_flux + end_flux
        
        total_flux.append(flux)
        
        
    return total_flux
    
       
       
 
def calcR_HK(H_flux,K_flux,V_flux,R_flux, b_v):
    
    '''Function to calculate the log(R_HK) value. The equations are taken from a paper by
    Mascareno et. al (2015) titled 'Rotation periods of late-type dwarf stars from time-series
    high-resolution spectroscopy of chromospheric indicators'
    '''
    
    #Calculate the S index
    
    alpha = 2.3*8.0
    
    S = alpha*(H_flux+K_flux)/(R_flux+V_flux)
    print S
    
    #Calculate C_cf
    
    power = 0.668- 1.270*(b_v) + 0.645*(b_v)**2 - 0.443*(b_v)**3

    C_cf = 10**power
    print C_cf
    
    #Calculate R_phot'
    
    power = 1.48*10**(-4)*np.exp(-4.3658*b_v)
    
    #R_phot = 10**power
    R_phot = power #Correction in the paper's formula
    print R_phot
    
    #Calculate log(R_HK)
    
    R_HK = (1.34*10**(-4)*C_cf*S) - R_phot
    
    return (R_HK)
    
        

    
    
     
    