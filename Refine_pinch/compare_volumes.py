from TricubicInterpolation import cTricubic as ti

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('..')

import kostas_filemanager as kfm
plt.style.use('kostas')

pinch1 = 'eclouds/refined_pinch1_cut_MTI1.0_MLI1.0_DTO1.0_DLO1.0.h5'
pinch2 = 'eclouds/refined_pinch1_cut_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5'
#pinch3 = 'eclouds/refined_Pinch7_cut_MTI2.0_MLI1.0_DTO1.0_DLO1.0.h5'
#epinch = 'eclouds/e_refined_Pinch7_cut_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5'

stats1 = kfm.h5_to_dict(pinch1, group='stats')
stats2 = kfm.h5_to_dict(pinch2, group='stats')
#stats3 = kfm.h5_to_dict(pinch3, group='stats')

bins1 = stats1['bins']
bins2 = stats2['bins']
#bins3 = stats3['bins']

fig = plt.figure(3,figsize=[18,5])
ax1 = fig.add_subplot(131)
ax1.hist(bins1[:-1], bins1, weights = stats1['log10_vx_hist'], histtype='step', log=True, color='r')
ax1.hist(bins2[:-1], bins2, weights = stats2['log10_vx_hist'], histtype='step', log=True, color='b')
#ax1.hist(bins3[:-1], bins3, weights = stats3['log10_vx_hist'], histtype='step', log=True, color='g')
ax1.set_xlabel('$\mathbf{\log_{10}\left(V_x\\right)}$')

ax2 = fig.add_subplot(132)
ax2.hist(bins1[:-1], bins1, weights = stats1['log10_vy_hist'], histtype='step', log=True, color='r')
ax2.hist(bins2[:-1], bins2, weights = stats2['log10_vy_hist'], histtype='step', log=True, color='b')
#ax2.hist(bins3[:-1], bins3, weights = stats3['log10_vy_hist'], histtype='step', log=True, color='g')
ax2.set_xlabel('$\mathbf{\log_{10}\left(V_y\\right)}$')

ax3 = fig.add_subplot(133)
ax3.hist(bins1[:-1], bins1, weights = stats1['log10_vz_hist'], histtype='step', log=True, color='r')
ax3.hist(bins2[:-1], bins2, weights = stats2['log10_vz_hist'], histtype='step', log=True, color='b')
#ax3.hist(bins3[:-1], bins3, weights = stats3['log10_vz_hist'], histtype='step', log=True, color='g')
ax3.set_xlabel('$\mathbf{\log_{10}\left(V_\zeta \\right)}$')

ax1.set_xlim(-6,0.5)
ax2.set_xlim(-6,0.5)
ax3.set_xlim(-6,0.5)
plt.show()




