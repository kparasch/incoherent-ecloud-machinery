import myfilemanager_sixtracklib as mfm
import numpy as np
import pickle
import matplotlib.pyplot as plt

in_folder = '/eos/user/k/kparasch/Decohered_STL_data/'
out_folder = '/eos/user/k/kparasch/processed_Decohered_STL_data/'
filelist10 = ['losses_sixtracklib_4271153.0.h5']
for ii in range(2,12):
     filelist10.append(f'losses_sixtracklib_4271153.{ii}.h5')
filelist20 = []
for ii in range(0,6):
    filelist20.append(f'losses_sixtracklib_4295944.{ii}.h5')
for ii in range(0,11):
    filelist20.append(f'losses_sixtracklib_4328392.{ii}.h5')
for ii in range(0,12):
    filelist20.append(f'losses_sixtracklib_4591293.{ii}.h5')
filelist15 = []
filelist15.append(f'losses_sixtracklib_4351523.0.h5')
filelist15.append(f'losses_sixtracklib_4351523.1.h5')
filelist15.append(f'losses_sixtracklib_4351523.2.h5')
filelist15.append(f'losses_sixtracklib_4351523.4.h5')
filelist15.append(f'losses_sixtracklib_4351523.5.h5')
for ii in range(0,12):
    filelist15.append(f'losses_sixtracklib_4395210.{ii}.h5')
filelist05 = []
for ii in range(0,12):
    filelist05.append(f'losses_sixtracklib_4448724.{ii}.h5')
filelist00 = []
for ii in range(0,8):
    filelist00.append(f'losses_sixtracklib_4487357.{ii}.h5')
for ii in range(9,12):
    filelist00.append(f'losses_sixtracklib_4487357.{ii}.h5')


filelists = [filelist00, filelist05, filelist10, filelist15, filelist20] 
fnames = ['collision_ptau0.0e-4', 'collision_ptau0.5e-4', 'collision_ptau1.0e-4', 'collision_ptau1.5e-4', 'collision_ptau2.0e-4'] 

def get_ob_input(myfile):
    print(myfile)
    #out_data       = mfm.h5_to_dict(myfile, group = 'output')
    in_data        = mfm.h5_to_dict(myfile, group = 'input')
    optics         = mfm.h5_to_dict(myfile, group = 'beam-optics')
    particle_on_CO = mfm.h5_to_dict(myfile, group = 'closed-orbit')
       
    
    e1 = optics['epsn_x'] / (optics['beta0'] * optics['gamma0'])
    e2 = optics['epsn_y'] / (optics['beta0'] * optics['gamma0'])
    
    invW = optics['invW']

    part_on_CO = np.array([particle_on_CO['x'],
                           particle_on_CO['px'],
                           particle_on_CO['y'],
                           particle_on_CO['py'],
                           particle_on_CO['zeta'],
                           particle_on_CO['delta']]).T

    X_init    = np.array( [in_data['init_x'].T,
                           in_data['init_px'].T,
                           in_data['init_y'].T,
                           in_data['init_py'].T,
                           in_data['init_zeta'].T,
                           in_data['init_delta'].T]).T
    
    X_init = np.array([X_init])
    X_init_norm = np.tensordot( X_init, invW, axes = (-1,1) )

    J1_init = 0.5*( X_init_norm[:,:,0]**2 + X_init_norm[:,:,1]**2)/e1
    J2_init = 0.5*( X_init_norm[:,:,2]**2 + X_init_norm[:,:,3]**2)/e2
    phi1_init = np.arctan2(-X_init_norm[:,:,1], X_init_norm[:,:,0])
    phi2_init = np.arctan2(-X_init_norm[:,:,3], X_init_norm[:,:,2])
    zeta_init = X_init[:,:,4]
    turn_init = np.zeros_like(J1_init)
    
    ob = { 'J1'   : J1_init,
           'J2'   : J2_init,
           'phi1' : phi1_init,
           'phi2' : phi2_init,
           'tau'  : zeta_init,
           'turn' : turn_init
         }

    return ob

def get_ob_first(myfile):
    print(myfile)
    out_data       = mfm.h5_to_dict(myfile, group = 'output')
    optics         = mfm.h5_to_dict(myfile, group = 'beam-optics')
    particle_on_CO = mfm.h5_to_dict(myfile, group = 'closed-orbit')
       
    
    e1 = optics['epsn_x'] / (optics['beta0'] * optics['gamma0'])
    e2 = optics['epsn_y'] / (optics['beta0'] * optics['gamma0'])
    
    invW = optics['invW']

    part_on_CO = np.array([particle_on_CO['x'],
                           particle_on_CO['px'],
                           particle_on_CO['y'],
                           particle_on_CO['py'],
                           particle_on_CO['zeta'],
                           particle_on_CO['delta']]).T

    X_first    = np.array([out_data['x_tbt_first'].T,
                           out_data['px_tbt_first'].T,
                           out_data['y_tbt_first'].T,
                           out_data['py_tbt_first'].T,
                           out_data['zeta_tbt_first'].T,
                           out_data['delta_tbt_first'].T]).T
    

    X_first_norm = np.tensordot( X_first - part_on_CO, invW, axes = (-1,1) )

    J1_first  = 0.5*( X_first_norm[:,:,0]**2 + X_first_norm[:,:,1]**2)/e1
    J2_first  = 0.5*( X_first_norm[:,:,2]**2 + X_first_norm[:,:,3]**2)/e2
    phi1_first = np.arctan2(-X_first_norm[:,:,1], X_first_norm[:,:,0])
    phi2_first = np.arctan2(-X_first_norm[:,:,3], X_first_norm[:,:,2])
    zeta_first = X_first[:,:,4]

    turn_first = out_data['at_turn_tbt_first']
    
    ob = { 'J1'   : J1_first,
           'J2'   : J2_first,
           'phi1' : phi1_first,
           'phi2' : phi2_first,
           'tau'  : zeta_first,
           'turn' : turn_first
         }

    return ob

def get_ob_skip(myfile):
    print(myfile)
    out_data       = mfm.h5_to_dict(myfile, group = 'output')
    optics         = mfm.h5_to_dict(myfile, group = 'beam-optics')
    particle_on_CO = mfm.h5_to_dict(myfile, group = 'closed-orbit')
       
    
    e1 = optics['epsn_x'] / (optics['beta0'] * optics['gamma0'])
    e2 = optics['epsn_y'] / (optics['beta0'] * optics['gamma0'])
    
    invW = optics['invW']

    part_on_CO = np.array([particle_on_CO['x'],
                           particle_on_CO['px'],
                           particle_on_CO['y'],
                           particle_on_CO['py'],
                           particle_on_CO['zeta'],
                           particle_on_CO['delta']]).T

    X_skip     = np.array([out_data['x_skip'].T,
                           out_data['px_skip'].T,
                           out_data['y_skip'].T,
                           out_data['py_skip'].T,
                           out_data['zeta_skip'].T,
                           out_data['delta_skip'].T]).T
    
    X_skip_norm = np.tensordot( X_skip - part_on_CO, invW, axes = (-1,1) )

    J1_skip = 0.5*( X_skip_norm[:,:,0]**2 + X_skip_norm[:,:,1]**2)/e1
    J2_skip = 0.5*( X_skip_norm[:,:,2]**2 + X_skip_norm[:,:,3]**2)/e2
    phi1_skip = np.arctan2(-X_skip_norm[:,:,1], X_skip_norm[:,:,0])
    phi2_skip = np.arctan2(-X_skip_norm[:,:,3], X_skip_norm[:,:,2])
    zeta_skip = X_skip[:,:,4]

    turn_skip = out_data['at_turn_skip']
    
    ob = { 'J1'   : J1_skip,
           'J2'   : J2_skip,
           'phi1' : phi1_skip,
           'phi2' : phi2_skip,
           'tau'  : zeta_skip,
           'turn' : turn_skip
         }

    return ob

def get_ob_last(myfile):
    print(myfile)
    out_data       = mfm.h5_to_dict(myfile, group = 'output')
    optics         = mfm.h5_to_dict(myfile, group = 'beam-optics')
    particle_on_CO = mfm.h5_to_dict(myfile, group = 'closed-orbit')
       
    
    e1 = optics['epsn_x'] / (optics['beta0'] * optics['gamma0'])
    e2 = optics['epsn_y'] / (optics['beta0'] * optics['gamma0'])
    
    invW = optics['invW']

    part_on_CO = np.array([particle_on_CO['x'],
                           particle_on_CO['px'],
                           particle_on_CO['y'],
                           particle_on_CO['py'],
                           particle_on_CO['zeta'],
                           particle_on_CO['delta']]).T

    X_last     = np.array([out_data['x_last'].T,
                           out_data['px_last'].T,
                           out_data['y_last'].T,
                           out_data['py_last'].T,
                           out_data['zeta_last'].T,
                           out_data['delta_last'].T]).T

    X_last_norm = np.tensordot( X_last - part_on_CO, invW, axes = (-1,1) )

    J1_last = 0.5*( X_last_norm[:,:,0]**2 + X_last_norm[:,:,1]**2)/e1
    J2_last = 0.5*( X_last_norm[:,:,2]**2 + X_last_norm[:,:,3]**2)/e2
    phi1_last = np.arctan2(-X_last_norm[:,:,1], X_last_norm[:,:,0])
    phi2_last = np.arctan2(-X_last_norm[:,:,3], X_last_norm[:,:,2])
    zeta_last = X_last[:,:,4]

    turn_last = out_data['at_turn_last']

    
    ob = { 'J1'   : J1_last,
           'J2'   : J2_last,
           'phi1' : phi1_last,
           'phi2' : phi2_last,
           'tau'  : zeta_last,
           'turn' : turn_last
         }

    return ob

def merge_obs(obs):
    ob_master = o
    
    J1 = obs[0]['J1']
    J2 = obs[0]['J2']
    phi1 = obs[0]['phi1']
    phi2 = obs[0]['phi2']
    tau = obs[0]['tau']
    turn = obs[0]['turn']
    for i in range(1,len(obs)):
        J1    = np.hstack((J1   , obs[i]['J1']))
        J2    = np.hstack((J2   , obs[i]['J2']))
        phi1  = np.hstack((phi1 , obs[i]['phi1']))
        phi2  = np.hstack((phi2 , obs[i]['phi2']))
        tau   = np.hstack((tau  , obs[i]['tau']))
        turn  = np.hstack((turn , obs[i]['turn']))

    ob = { 'J1'   : J1,
           'J2'   : J2,
           'phi1' : phi1,
           'phi2' : phi2,
           'tau'  : zeta,
           'turn' : turn
         }
    return ob


filelist = filelists[0]
suffix = '_first.h5'
fname = fnames[0]+suffix
print(fname)
outfile = out_folder + fname
obs = []
for myfile in filelist:
    obs.append( get_ob_first(in_folder+myfile) )
    master_ob = merge_obs(obs)
optics         = kfm.h5_to_dict(myfile, group = 'beam-optics')
particle_on_CO = kfm.h5_to_dict(myfile, group = 'closed-orbit')
kfm.dict_to_h5(master_ob, outfile, group='data', readwrite_opts='w')
kfm.dict_to_h5(optics, outfile, group='optics', readwrite_opts='a')
kfm.dict_to_h5(particle_on_CO, outfile, group='closed-orbit', readwrite_opts='a')



