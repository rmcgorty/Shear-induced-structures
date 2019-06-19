# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 14:02:33 2019

@author: rmcgorty
"""

import numpy as np
from numpy.fft import *
import glob
import matplotlib
from matplotlib import pyplot as plt
import pickle
import sys

#Add the following to path. This is for importing the 'tiff_file' module.
#  You can change this location to some other where 'tiff_file.py' is located. 
sys.path.append("C:\\Users\\RMCGORTY\\Documents\\GitHub\\Differential-Dynamic-Microscopy---Python")

import tiff_file

#Specify the data directory
data_dir = "Z:\\Michelle_research\\Data\\2019-06-13_shearing\\test_10\\"

#This will find all filenames matching the pattern given
tif_files = glob.glob(data_dir + "*.tif")

#Report out how many tiff files are found in that directory
print "Found %i tiff files." % len(tif_files)


#Loop through all images in the directory
for index, filename in enumerate(tif_files):
    image = tiff_file.imread(filename).sum(axis=-1) #Load the images (the 'sum' is for converting RGB to monochrome)
    
    #Next two lines just normalize the image data
    data = image - image.mean()
    data = data / data.std()
    
    #To compute the autocorrelation, use fft
    res = np.real(fftshift(ifft2(fft2(data)*np.conj(fft2(data)))))
    res = res / (data.shape[0]*data.shape[1])
    
    if index==0:
        res_all = res
    else:
        res_all += res

#Divide by total images so we have an average autocorrelation of all frames
res_all = res_all / len(tif_files)

#Now you can use 'plt.matshow(res_all)' to see the autocorrelation.
#Or 'plt.plot(res_all[400,400:],'-ro')' to a section along the x-axis. 