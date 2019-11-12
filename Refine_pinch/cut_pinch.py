import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

import sys
sys.path.append('..')
import kostas_filemanager as kfm

pinch_folder = '/home/kparasch/Storage/MyPinches/'
pinch = 'Pinch7'
print('Loading file...')
ob = kfm.h5_to_dict(pinch_folder + pinch +'.h5')

print('Shape of xg',ob['xg'].shape)
print('Shape of yg',ob['yg'].shape)
print('Shape of zg',ob['zg'].shape)
print('Shape of phi',ob['phi'].shape)

i_zero = np.argmin(np.abs(ob['xg']))
j_zero = np.argmin(np.abs(ob['yg']))

Nx = ob['phi'].shape[0]
Ny = ob['phi'].shape[1]
Nz = ob['phi'].shape[2]

N_keep_x = 101
N_keep_y = 101
N_keep_z = 200

x0 = (Nx - N_keep_x)//2
y0 = (Ny - N_keep_y)//2
z0 = (Nz - N_keep_z)//2

dict_new_file = {}

dict_new_file['phi'] = ob['phi'][x0:-x0, y0:-y0, z0:-z0]
dict_new_file['rho'] = ob['rho'][x0:-x0, y0:-y0, z0:-z0]
dict_new_file['xg' ] = ob['xg'][x0:-x0]
dict_new_file['yg' ] = ob['yg'][y0:-y0]
dict_new_file['zg' ] = ob['zg'][z0:-z0]

print('Old shape: ', ob['phi'].shape)
print('New shape: ', dict_new_file['phi'].shape)
print('x, y, z new Ranges:')
print(dict_new_file['xg'][0], dict_new_file['xg'][-1])
print(dict_new_file['yg'][0], dict_new_file['yg'][-1])
print(dict_new_file['zg'][0], dict_new_file['zg'][-1])

kfm.dict_to_h5(dict_new_file, pinch_folder + pinch + '_cut.h5', compression_opts=0)
