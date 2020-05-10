import numpy as np
import shutil

import sys
sys.path.append('../Tools')
import kostas_filemanager as kfm

list_of_pinches = ['/eos/user/k/kparasch/Pinches/ArcDipole_07.h5',
                   '/eos/user/k/kparasch/Pinches/ArcDipole_07_2000.h5'
                  ]
numbers_of_pinches = [4000,2000]

temp_file = 'temp.h5'
out_file = '/eos/user/e/ecincohe/Raw_Pinches/LHC_ArcDip_1.35sey_0.7e11ppb.h5'

fname0 = list_of_pinches[0]
grid_dict = kfm.h5_to_dict(fname0, group='grid')
stats_dict = kfm.h5_to_dict(fname0, group='stats')
kfm.dict_to_h5(grid_dict, temp_file, group='grid', readwrite_opts='w')
kfm.dict_to_h5(stats_dict, temp_file, group='stats', readwrite_opts='a')

nx = len(grid_dict['xg'])
ny = len(grid_dict['yg'])

result_dict = {'phi': np.zeros([nx,ny]), 'rho': np.zeros([nx,ny])}

n_slices = len(grid_dict['zg'])

total_pinches = sum(numbers_of_pinches)

for kk in range(n_slices):
    print(f'Slice: {kk}/{n_slices}')
    this_slice = f'slices/slice{kk}'
    phi = np.zeros([nx,ny])
    result_dict['phi'][:,:] = 0.
    result_dict['rho'][:,:] = 0.
    for fname, NN in zip(list_of_pinches, numbers_of_pinches):
        this_dict = kfm.h5_to_dict(fname, group=this_slice)
        result_dict['phi'] += this_dict['phi']*((1.*NN)/(1.*total_pinches))
        result_dict['rho'] += this_dict['rho']*((1.*NN)/(1.*total_pinches))
    kfm.dict_to_h5(result_dict, temp_file, group=this_slice, readwrite_opts='a')


shutil.copy(temp_file, out_file)

