import sys
import shutil
#import contextlib2 as contextlib #pip install contextlib2 in python2 or use contextlib in python3.4>
import contextlib
import os
import numpy as np
import time
import multiprocessing

#sys.path.append(os.path.abspath('../Tools'))
#sys.path.append(os.path.abspath('../../packages/PyFRIENDS'))
import kostas_filemanager as kfm
import PyPARIS_sim_class.Simulation as sim_mod
import Simulation_parameters as pp

pinch_name = sys.argv[1]
workers = int(sys.argv[2])
n_pinches_to_average = int(sys.argv[3])
out_dir = sys.argv[4]
#workers=2
#n_pinches_to_average = 4
#n_pinches_to_average = 4000
#save_efields = True
#save_rho = True
#transpose_to_natural_ordering_of_xyz = True

pp_dict = {}
for attr in dir(pp):
    if attr[0] != '_':
        temp = getattr(pp, attr)
        if temp is None:
            temp = 'None'
        pp_dict[attr] = temp
pp_dict['pinch_name'] = pinch_name
pp_dict['max_workers'] = workers
pp_dict['n_pinches_to_average'] = n_pinches_to_average
pp_dict['out_dir'] = out_dir


def one_pinch(mydict, lock, N_pinches=1, save_sigmas_and_coords=False, idd=0, grid=None):
    os.mkdir('temp'+str(idd))
    shutil.copytree('pyecloud_config', 'temp'+str(idd)+'/pyecloud_config')
    shutil.copyfile('Simulation_parameters.py', 'temp'+str(idd)+'/Simulation_parameters.py')
    shutil.copyfile('LHC_chm_ver.mat', 'temp'+str(idd)+'/LHC_chm_ver.mat')
    os.chdir('temp'+str(idd))
    if os.path.exists('simulation_status.sta'):
        os.remove('simulation_status.sta')
    
    ring = sim_mod.get_serial_CPUring()
    
    list_slice_objects = ring.pieces_to_be_treated # Head is the last element
    list_machine_elements = ring.sim_content.mypart
    # Truck the beam to the first ecloud (excluded)
    for ee in list_machine_elements:
        if ee in ring.sim_content.my_list_eclouds:
            first_ecloud = ee
            break
        for ss in list_slice_objects[::-1]:
            ee.track(ss)
    
    # Record pinch info
    first_ecloud.save_ele_distributions_last_track = True
#    first_ecloud.save_ele_field = True
    first_ecloud.save_ele_potential = True
    
    first_ecloud._reinitialize() # Needed to prepare storage space
    
    N_slices = len(list_slice_objects)
    z_centers = []
    t_start = time.mktime(time.localtime())
    for i_ss, ss in enumerate(list_slice_objects[::-1]):
        if np.mod(i_ss, 20)==0:
            print("%d / %d"%(i_ss, N_slices))
        first_ecloud.track(ss)
        if ss.slice_info != 'unsliced':
            z_centers.append(ss.slice_info['z_bin_center'])
    
    first_ecloud._finalize()
    z_centers = z_centers[::-1] # HEADTAIL convention
    
    t_end = time.mktime(time.localtime())
    
    print('Track time %.2f s' % (t_end - t_start))

    while 1:
        if len([temp_dirs for temp_dirs in os.listdir('.') if 'temp' in temp_dirs]) < 110:
            break
        else:
            time.sleep(60)
   
    if save_sigmas_and_coords:
        grid['sigma_x_beam'] = ring.sim_content.bunch.sigma_x()
        grid['sigma_y_beam'] = ring.sim_content.bunch.sigma_y()
        grid['sigma_z_beam'] = ring.sim_content.bunch.sigma_z()
        grid['xg'] = first_ecloud.spacech_ele.xg
        grid['yg'] = first_ecloud.spacech_ele.yg
        grid['zg'] = z_centers

    first_ecloud.phi_ele_last_track /= (1.*N_pinches)
    first_ecloud.rho_ele_last_track /= (1.*N_pinches)
    lock.acquire()
    if 'phi' in mydict.keys():
        mydict['phi'] += first_ecloud.phi_ele_last_track
        mydict['rho'] += first_ecloud.rho_ele_last_track
    else:
        mydict['phi'] = first_ecloud.phi_ele_last_track
        mydict['rho'] = first_ecloud.rho_ele_last_track
    lock.release()

    os.chdir('..')
  
    return idd

def kern(i, mydict, lock):
    print(f'Running: {i+1}')
    with open('stdout'+str(i)+'.out','w') as f:
        with contextlib.redirect_stdout(f):
            return one_pinch(mydict, lock, N_pinches=n_pinches_to_average, save_sigmas_and_coords=False, idd=i+1, grid=None)

manager = multiprocessing.Manager()
result_dict = manager.dict()
grid_dict = manager.dict()
lock = manager.Lock()

i0 = one_pinch(result_dict, lock, N_pinches=n_pinches_to_average, save_sigmas_and_coords=True, grid=grid_dict)

with multiprocessing.Pool(workers) as pool:
    result_list = pool.starmap_async(kern, [(ii, result_dict, lock) for ii in range(n_pinches_to_average-1)])
    print(result_list.get())


# Save pinch to file

out_pinch = pinch_name + '.h5'

kfm.dict_to_h5(grid_dict, out_pinch, group='grid', readwrite_opts='w')
for i in range(result_dict['phi'].shape[0]):
    kfm.dict_to_h5({'phi' : result_dict['phi'][i,:,:], 'rho' : result_dict['rho'][i,:,:]}, out_pinch, group='slices/slice%d'%i, readwrite_opts='a')

kfm.dict_to_h5({'ti_method' : 'FD'}, out_pinch, group='stats', readwrite_opts='a')
kfm.dict_to_h5({'sim' : str(pp_dict)}, out_pinch, group='Simulation_parameters', readwrite_opts='a')

final_destination = out_dir+'/'+pinch_name+'.h5'
shutil.copyfile(out_pinch, final_destination)
print(f"max x : {grid_dict['xg'][-1]/grid_dict['sigma_x']}sigmas")
print(f"max y : {grid_dict['yg'][-1]/grid_dict['sigma_y']}sigmas")
print(f"max z : {grid_dict['zg'][-1]/grid_dict['sigma_z']}sigmas")
print(f"Sigma_x: {grid_dict['sigma_x']}")
print(f"Sigma_y: {grid_dict['sigma_y']}")
print(f"Sigma_z: {grid_dict['sigma_z']}")
print(f"Pinch can be found in {final_destination}")
