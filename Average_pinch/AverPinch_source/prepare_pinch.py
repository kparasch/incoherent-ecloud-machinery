import sys
import shutil
#import contextlib2 as contextlib #pip install contextlib2 in python2 or use contextlib in python3.4>
import contextlib
import os
import numpy as np
import time
import concurrent.futures

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
save_rho = True
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

#os.mkdir(pinch_name)
#shutil.copyfile('AP_source/Simulation_parameters.py', pinch_name+'/Simulation_parameters.py')
#shutil.copytree('AP_source/Simulation_parameters.py', pinch_name+'/')
#sys.path.append(os.path.abspath(pinch_name))
#os.chdir(pinch_name)
#import Simulation_parameters as pp

executor = concurrent.futures.ProcessPoolExecutor(max_workers=workers)




def one_pinch(save_rho=True, save_sigmas=False, save_coords=False,idd=0):
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

    dd = {}   
    dd['phi'] = first_ecloud.phi_ele_last_track
    if save_rho:
        dd['rho'] = first_ecloud.rho_ele_last_track
    if save_sigmas:
        dd['sigma_x_beam'] = ring.sim_content.bunch.sigma_x()
        dd['sigma_y_beam'] = ring.sim_content.bunch.sigma_y()
        dd['sigma_z_beam'] = ring.sim_content.bunch.sigma_z()
    if save_coords:
        dd['xg'] = first_ecloud.spacech_ele.xg
        dd['yg'] = first_ecloud.spacech_ele.yg
        dd['zg'] = z_centers

    kfm.dict_to_h5(dd,'temp_pinch.h5')
    os.chdir('..')
    #os.rmdir('temp'+str(idd))
  
    return idd

def kern(i):
    with open('stdout'+str(i)+'.out','w') as f:
        with contextlib.redirect_stdout(f):
            return one_pinch(save_rho=save_rho, save_sigmas=False, save_coords=False,idd=i+1)

i0 = one_pinch(save_rho=save_rho, save_sigmas=True, save_coords=True)
dd = kfm.h5_to_dict('temp'+str(i0)+'/temp_pinch.h5')
os.remove('temp'+str(i0)+'/temp_pinch.h5')

grid_dict = { 'xg' : dd['xg'],
              'yg' : dd['yg'],
              'zg' : dd['zg'],
              'sigma_x' : dd['sigma_x_beam'], 
              'sigma_y' : dd['sigma_y_beam'], 
              'sigma_z' : dd['sigma_z_beam']
            }

dd['phi'] /= 1.*n_pinches_to_average
if save_rho:
    dd['rho'] /= 1.*n_pinches_to_average

results=[]
for i in range(n_pinches_to_average-1):
    results.append(executor.submit(kern,i))

jj = 0
while len(results):
    for j in range(len(results)):
        if results[j].done():
            print('something is done')
            idd = results[j].result()
            dd_temp = kfm.h5_to_dict('temp'+str(idd)+'/temp_pinch.h5')
            jj += 1
            dd['phi'] += dd_temp['phi']/(1.*n_pinches_to_average)
            if save_rho:
                dd['rho'] += dd_temp['rho']/(1.*n_pinches_to_average)
            os.remove('temp'+str(idd)+'/temp_pinch.h5')
            if int(idd) > 10:
                os.remove('stdout'+str(idd)+'.out')
                shutil.rmtree('temp'+str(idd))
            del dd_temp
            del results[j]
            print('Finished #%d/%d pinches.'%(jj+1,n_pinches_to_average))
            break
       
#for i in range(n_pinches_to_average-1):
#    print('Finished #%d/%d pinches.'%(i,n_pinches_to_average))
#    #print('Running pinch #%d/%d.'%(i,n_pinches_to_average))
#    dd_temp = one_pinch(save_rho=save_rho, save_sigmas=False, save_coords=False)

#if save_efields:
#    dx = dd['xg'][1] - dd['xg'][0]
#    dy = dd['yg'][1] - dd['yg'][0]
#    dz = dd['zg'][1] - dd['zg'][0]
#    dd['Ex'] = np.zeros_like(dd['phi'])
#    dd['Ey'] = np.zeros_like(dd['phi'])
#    dd['Ez'] = np.zeros_like(dd['phi'])
#    
#    dd['Ex'][:,1:-1,:] = -0.5/dx*( dd['phi'][:,2:,:] - dd['phi'][:,0:-2,:])
#    dd['Ey'][:,:,1:-1] = -0.5/dy*( dd['phi'][:,:,2:] - dd['phi'][:,:,0:-2])
#    dd['Ez'][1:-1,:,:] = -0.5/dz*( dd['phi'][2:,::] - dd['phi'][0:-2,:,:])

#if transpose_to_natural_ordering_of_xyz:
#    dd['phi'] = dd['phi'].transpose(1,2,0)
#    dd['rho'] = dd['rho'].transpose(1,2,0)
#    dd['Ex']  = dd['Ex'].transpose(1,2,0)
#    dd['Ey']  = dd['Ey'].transpose(1,2,0)
#    dd['Ez']  = dd['Ez'].transpose(1,2,0)

# Save pinch to file

final_destination = out_dir+'/'+pinch_name+'.h5'

kfm.dict_to_h5(grid_dict, final_destination, group='grid', readwrite_opts='w')
for i in range(dd['phi'].shape[0]):
    kfm.dict_to_h5({'phi' : dd['phi'][i,:,:], 'rho' : dd['rho'][i,:,:]}, final_destination, group='slices/slice%d'%i, readwrite_opts='a')

kfm.dict_to_h5({'ti_method' : 'FD'}, final_destination, group='stats', readwrite_opts='a')
kfm.dict_to_h5(pp_dict, final_destination, group='Simulation_parameters', readwrite_opts='a')
print(f"max x : {grid_dict['xg'][-1]/grid_dict['sigma_x']}sigmas")
print(f"max y : {grid_dict['yg'][-1]/grid_dict['sigma_y']}sigmas")
print(f"max z : {grid_dict['zg'][-1]/grid_dict['sigma_z']}sigmas")
print(f"Sigma_x: {grid_dict['sigma_x']}")
print(f"Sigma_y: {grid_dict['sigma_y']}")
print(f"Sigma_z: {grid_dict['sigma_z']}")
print(f"Pinch can be found in {final_destination}")
