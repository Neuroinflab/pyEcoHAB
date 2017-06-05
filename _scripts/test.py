# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:32:16 2016

@author: Janek
"""
import numpy as np
import pylab as plt
import scipy.stats as st
dane = np.load('result.npy')
error= np.load('result_err.npy')
maxi = np.max(dane)
min_err = st.scoreatpercentile(error, 2.5, axis=2)
max_err = st.scoreatpercentile(error, 97.5, axis=2)
plt.imshow(min_err,cmap=plt.gray(),interpolation='none',vmin=0, vmax=maxi)
plt.show()
plt.imshow(dane[:,:,0],cmap=plt.gray(),interpolation='none',vmin=0, vmax=maxi)
plt.show()
plt.imshow(max_err,cmap=plt.gray(),interpolation='none',vmin=0, vmax=maxi)
plt.show()

plt.hist(error.flatten(),range=[5, 25],bins= 20)
plt.show()
#print error[3,8,:]