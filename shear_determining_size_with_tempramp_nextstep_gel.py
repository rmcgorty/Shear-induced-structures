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
from matplotlib import cm
import pickle
import sys

font_plt = {'family': 'serif',
            'color':  'k',
            'weight': 'normal',
            'size': 10,}

txt_plt = {'family': 'serif',
            'color':  'k',
            'weight': 'normal',
            'size': 7,}


# useful colors
gray10 = '#191919'
gray20 = '#333333'
gray30 = '#4C4C4C'
gray40 = '#666666'
gray50 = '#7F7F7F'
gray60 = '#999999'
gray70 = '#B2B2B2'
gray80 = '#CCCCCC'
gray90 = '#E5E5E5'

#Add the following to path. This is for importing the 'tiff_file' module.
#  You can change this location to some other where 'tiff_file.py' is located. 
sys.path.append("C:\\Users\\RMCGORTY\\Documents\\GitHub\\Differential-Dynamic-Microscopy---Python")

import tiff_file

#FOR SHEARING SCOPE: 0.497 microns per pixel

#Specify the data directory and data file
data_dir = "Z:\\Michelle_research\\Data\\2019-06-17_shearing\\23_osc_100gel_cooling_2hz\\"
data_file = "23_osc_100gel_cooling_rollingbgsubtract-200.tif"

#data = pickle.load(open(data_dir+data_file[:-4]+"_Analysis.p",'rb'))

results = data['corr_results']
frames = data['frames_used']

fig = plt.figure(figsize=(8,6), constrained_layout=True)

linew = 7.0

#frames_to_plot = [5,7,9,11,13]
frames_to_plot = [4,7,10,13]
colors = [cm.cubehelix(x) for x in np.linspace(0.05, 0.7, len(frames_to_plot)) ]
colors = [cm.gray(x) for x in np.linspace(0.05, 0.65, len(frames_to_plot)) ]

x = np.arange(0,400)*0.497
for i in range(len(frames_to_plot)):
    plt.plot(x[::1],results[frames_to_plot[i],400,400::1],linestyle='-',marker='.',c=colors[i],lw=linew)

plt.xlabel('Distance (microns)')
plt.ylabel('Autocorrelation')


#Specify the data directory and data file
data_dir = "Z:\\Michelle_research\Data\\2019-06-15_shearing\\shear_0.6_cooling\\"
data_file = "shear_0.6_cooling_rollingballbgsubtract-50.tif"

#data_c = pickle.load(open(data_dir+data_file[:-4]+"_Analysis.p",'rb'))

results_c = data_c['corr_results']
frames_c = data_c['frames_used']

fig = plt.figure(figsize=(2.75,2.25), constrained_layout=True)

linew = 2.5



frames_to_plot = [4,6,8,10,12]
colors = [cm.gist_heat(x) for x in np.linspace(0.05, 0.9, len(frames_to_plot)) ]
x = np.arange(0,400)*0.497
for i in range(len(frames_to_plot)):
    plt.plot(x[::1],1+results[frames_to_plot[i],400,400::1],linestyle='-',marker='None',c=colors[i],lw=linew)

frames_to_plot = [20,27,34,41]
colors = [cm.gist_heat(x) for x in np.linspace(0.05, 0.9, len(frames_to_plot)) ]
colors.reverse()
x = np.arange(0,400)*0.497
for i in range(len(frames_to_plot)):
    plt.plot(x[::1],results_c[frames_to_plot[i],400,400::1],linestyle='-',marker='None',c=colors[i],lw=linew)

plt.xlim(0,7)
plt.xlabel('Distance (microns)')
plt.ylabel('Autocorrelation')

plt.savefig('fig_forproposal.svg', dpi=600)


plt.figure(figsize=(2,2.25),constrained_layout=True)
plt.matshow(results[:,390:411,400:500].mean(axis=-2)[0:25,0:30],fignum=0)
plt.savefig(data_dir+'heating_0.6.svg', dpi=600)

plt.figure(figsize=(2,2.25),constrained_layout=True)
plt.matshow(results_c[:,390:411,400:500].mean(axis=-2)[-25:,0:30],fignum=0)
plt.savefig(data_dir+'cooling_0.6.svg', dpi=600)

'''
fig = plt.figure(figsize=(8,6), constrained_layout=True)

linew = 7.0

#frames_to_plot = [5,7,9,11,13]
frames_to_plot = [20,25,30,35,40,45]
colors = [cm.cubehelix(x) for x in np.linspace(0.05, 0.7, len(frames_to_plot)) ]
colors = [cm.gray(x) for x in np.linspace(0.05, 0.65, len(frames_to_plot)) ]

x = np.arange(0,400)*0.497
for i in range(len(frames_to_plot)):
    plt.plot(x[::1],results_c[frames_to_plot[i],400,400::1],linestyle='-',marker='.',c=colors[i],lw=linew)

plt.xlabel('Distance (microns)')
plt.ylabel('Autocorrelation')
'''

#allResults = {}
#allResults['std_threshold'] = std_threshold
#allResults['num_to_avg'] = num_to_avg
#allResults['data_file'] = data_file
#allResults['frames_used'] = np.array(frames_used)
#allResults['data_dir'] = data_dir
#allResults['fft_results'] = fft_res
#allResults['corr_results'] = results
#pickle.dump(allResults, open(data_dir+data_file[:-4]+"_Analysis.p",'wb'))



################################################################################
# OLD STUFF
################################################################################
'''
#Specify the data directory and data file
data_dir = "Z:\\Michelle_research\Data\\2019-06-14_shearing\\shear_2.5_cooling\\"
data_file = "shear_2.5_cooling.tif"

image = tiff_file.imread(data_dir+data_file)[:,:256,:256]
print "Image has %i frames." % image.shape[0]

stds = image.std(axis=-1).std(axis=-1) #find the std deviation of each image
#This stds should be a vector of length equal to the number of frames

std_threshold = 0.5 #PLOT THE STD DEVS TO SEE WHERE THRESHOLD SHOULD BE!
w = np.where(stds<std_threshold)

num_time_pts = 30
results = np.zeros((num_time_pts,image.shape[1],image.shape[2]))

j = 0 #initialize counter which will increment during for loop

frames_used = []


for i in np.linspace(0,image[w[0]].shape[0],num_time_pts,endpoint=False,dtype=int):
    frames_used.append(w[0][i])
    for k in range(3):

        #Next two lines just normalize the image data
        data = image[w[0]][i+k] - image[w[0]][i+k].mean()
        data = data / data.std()
        
        #To compute the autocorrelation, use fft
        res = np.real(fftshift(ifft2(fft2(data)*np.conj(fft2(data)))))
        res = res / (data.shape[0]*data.shape[1])
        
        results[j] += res
    results[j] = results[j]/100.0
    j+=1

'''
