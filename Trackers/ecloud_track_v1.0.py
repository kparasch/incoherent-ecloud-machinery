import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
import shutil
import copy
import time
import sys

sys.path.append('../Tools')
import distribution
import ecloud_sixtracklib_helpers as ec_stl
import kostas_filemanager as kfm

import pysixtrack
import sixtracklib
import NAFFlib


parser = argparse.ArgumentParser(description='Tracking simulations with e-clouds')
parser.add_argument('--jobtype', nargs='?', default=None, type=str)
parser.add_argument('--device', nargs='?', default=None, type=str)
parser.add_argument('--copy_destination', nargs='?', default=None, type=str)
parser.add_argument('--ecloud_strength', nargs='?', default=1., type=float)
parser.add_argument('--simulation_input', nargs='?', default='simulation_input.pkl', type=str)
parser.add_argument('--eclouds_folder', nargs='?', default='', type=str)
parser.add_argument('--ptau_max', nargs='?', default=0., type=float)
parser.add_argument('--max_tau', nargs='?', default=0.1, type=float)
parser.add_argument('--output', nargs='?', default='temp.h5', type=str)
parser.add_argument('--skip_turns', nargs='?', default=10000, type=int)
parser.add_argument('--turns_per_checkpoint', nargs='?', default=1000000, type=int)
parser.add_argument('--DA_turns', nargs='?', default=1000000, type=int)
parser.add_argument('--DA_r_N', nargs='?', default=200, type=int)
parser.add_argument('--DA_theta_N', nargs='?', default=100, type=int)
parser.add_argument('--last_checkpoint', nargs='?', default=10, type=int)
parser.add_argument('--n_particles', nargs='?', default=20000, type=int)
parser.add_argument('--n_sigma', nargs='?', default=5.5, type=float)
parser.add_argument('--seed', nargs='?', default=0, type=int)
parser.add_argument('--mb', nargs='?', default=None, type=str)
parser.add_argument('--mqf', nargs='?', default=None, type=str)
parser.add_argument('--mqd', nargs='?', default=None, type=str)

args = parser.parse_args()

if args.jobtype not in ['FMA', 'DA', 'LE']:
    raise Exception(f'Job type unknown: {args.jobtype}')

if args.copy_destination is not None:
    if args.copy_destination[-1] != '/':
        args.copy_destination += '/'

args_dict = vars(args).copy()
for key in args_dict.keys():
    if args_dict[key] is None:
        args_dict[key] = 'None'

allowed_eclouds = ['mb', 'mqf', 'mqd']
ecloud_sources = {key:args.eclouds_folder+getattr(args,key) for key in allowed_eclouds if getattr(args,key)}
print('E-clouds enabled:')
for key in ecloud_sources.keys():
    print(f'{key}:{ecloud_sources[key]}')


with open(args.simulation_input, 'rb') as fid:
    sim_input = pickle.load(fid)
sim_input['line'] = pysixtrack.Line.from_dict(sim_input['line'], keepextra=True)
optics = sim_input['optics']

for key in sim_input['eclouds_info']['length'].keys():
    sim_input['eclouds_info']['length'][key] *= args.ecloud_strength/(optics['beta0']*optics['p0c_eV'])
        
if ecloud_sources:
    old_optics = copy.deepcopy(optics)
    new_optics = ec_stl.update_optics(sim_input, ecloud_sources)
    sim_input['optics'] = new_optics
    optics = new_optics

start_time = time.time()


if args.jobtype == 'LE':
    if args.copy_destination is not None:
        shutil.copy(args.copy_destination + args.output, args.output)

    checkpoint = kfm.h5_to_dict(args.output, group='checkpoint')['checkpoint']
    print(f'Starting from checkpoint: {checkpoint}!')

    n_stores1 = int(args.turns_per_checkpoint/args.skip_turns)
    n_stores2 = 1
    epsn_1 = 2.0e-6
    epsn_2 = 2.0e-6
    
    sim_input['line'].append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores1, start=checkpoint*args.turns_per_checkpoint + args.skip_turns-1, skip=args.skip_turns, is_rolling=True),'monitor1')
    sim_input['line'].append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores2, start=checkpoint*args.turns_per_checkpoint, is_rolling=True),'monitor2')

    if checkpoint == 0:
        init_denormalized_6D = distribution.get6D_with_matched_J3_shell(n_particles=args.n_particles, 
                                                                        n_sigma=args.n_sigma, 
                                                                        ptau_max=args.ptau_max, 
                                                                        epsn_1=epsn_1, 
                                                                        epsn_2=epsn_2,
                                                                        optics=optics, 
                                                                        seed=args.seed
                                                                       )
    
        init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, sim_input['partCO'])
    
        init_dict = {'x': init_denormalized_6D[:,0], 'px': init_denormalized_6D[:,1],
                     'y': init_denormalized_6D[:,2], 'py': init_denormalized_6D[:,3],
                     'tau': init_denormalized_6D[:,4], 'ptau': init_denormalized_6D[:,5]
                    }
    
        ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])
    else:
        ps = distribution.load_state(args.output, p0c_eV=optics['p0c_eV'])

elif args.jobtype == 'DA':
    n_stores = 1
    turn_to_track = args.DA_turns
    n_particles_approx = 20000
    #n_sigma = 5.7
    epsn_1 = 3.5e-6 #
    epsn_2 = 3.5e-6
    
    se1 = np.sqrt(epsn_1/optics['gamma0']/optics['beta0'])
    se2 = np.sqrt(epsn_2/optics['gamma0']/optics['beta0'])
    
    sim_input['line'].append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores,is_rolling=True),'monitor1')
    init_denormalized_6D, A1_A2_in_sigma, da_particles = distribution.get_DA_distribution(
                                                         n_sigma=args.n_sigma, ptau_max=args.ptau_max, 
                                                         epsn_1=epsn_1, epsn_2=epsn_2, 
                                                         optics=optics,
                                                         r_N=args.DA_r_N, theta_N=args.DA_theta_N
                                                         
                                                                                         )
    
    init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, sim_input['partCO'])
    ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])
    
    J1, J2 = distribution.J1_J2_from_physical(init_denormalized_6D, optics['invW'], sim_input['partCO'])

elif args.jobtype == 'FMA':
    n_stores = 1000
    #n_sigma = 0.001#5.6
    epsn_1 = 3.5e-6
    epsn_2 = 3.5e-6
    
    se1 = np.sqrt(epsn_1/optics['gamma0']/optics['beta0'])
    se2 = np.sqrt(epsn_2/optics['gamma0']/optics['beta0'])
    
    sim_input['line'].append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores,is_rolling=True),'monitor1')
    
    init_denormalized_6D = distribution.get_fma_distribution(n_particles=args.n_particles, 
                                                             n_sigma=args.n_sigma, ptau_max=args.ptau_max, 
                                                             epsn_1=epsn_1, epsn_2=epsn_2, 
                                                             optics=optics, seed=args.seed
                                                            )
    
    init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, sim_input['partCO'])
    
    ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])
    
    J1, J2 = distribution.J1_J2_from_physical(init_denormalized_6D, optics['invW'], sim_input['partCO'])


ecloud_lattice = ec_stl.LatticeWithEclouds(sim_input, ps, ecloud_types=list(ecloud_sources.keys()), device=args.device)

for key in ecloud_sources.keys():
    ecloud_lattice.add_tricub_data(ecloud_sources[key], key, max_z=args.max_tau)
ecloud_lattice.remove_dipolar_kicks()
ecloud_lattice.finalize_assignments()
job = ecloud_lattice.job

end_setup_time = time.time()
print(f'Setting up time: {(end_setup_time - start_time)/60.:.4f}mins')

if args.jobtype =='LE':
    start_tracking = time.time()
    for cc in range(checkpoint+1, args.last_checkpoint+1):
        checkpoint_dict = {'checkpoint': cc} 
        skip_dicts, last_dict = ec_stl.track_to_checkpoint(job, 
                              n_particles=args.n_particles, checkpoint=cc, 
                              checkpoint_turns=args.turns_per_checkpoint, 
                              monitor1_stores=n_stores1, monitor2_stores=n_stores2, 
                              skip_turns=args.skip_turns
                                                          )
    
    
        if cc == 1:
            kfm.dict_to_h5(init_dict, args.output, group='init', readwrite_opts='a')
            kfm.dict_to_h5(optics, args.output, group='optics', readwrite_opts='a')
            kfm.dict_to_h5(sim_input['partCO'], args.output, group='closed-orbit', readwrite_opts='a')
            kfm.dict_to_h5(last_dict, args.output, group='last', readwrite_opts='a')
            time_dict = {'setup_time_mins' : 0., 'tracking_time_mins' : 0. }
            kfm.dict_to_h5(time_dict, args.output, group='time', readwrite_opts='a')
            kfm.dict_to_h5(args_dict, args.output, group='args', readwrite_opts='a')
        else:
            kfm.overwrite(last_dict, args.output, group='last')
    
        for (turn, skip_dict) in skip_dicts:
            kfm.dict_to_h5(skip_dict, args.output, group=f'turn{turn:d}', readwrite_opts='a')
    
        kfm.overwrite(checkpoint_dict, args.output, group='checkpoint')
        kfm.dict_to_h5(checkpoint_dict, args.output, group=f'checkpoint{cc}')
    
        if args.copy_destination is not None:
            shutil.copy(args.output, args.copy_destination + args.output)
    
        print(f'Reached checkpoint: {cc}/{args.last_checkpoint}')
    end_tracking = time.time()
    print(f'Tracking time: {(end_tracking - start_tracking)/60.:.4f}mins')
    
    time_dict = {'setup_time_mins' : (end_setup_time - start_time)/60.,
                 'tracking_time_mins' : (end_tracking - start_tracking)/60.
                }
    
    kfm.overwrite(args_dict, args.output, group='args')
    kfm.overwrite(time_dict, args.output, group='time')

    if args.copy_destination is not None:
        shutil.copy(args.output, args.copy_destination + args.output)

elif args.jobtype == 'DA':
    start_tracking = time.time()
    job.track_until(turn_to_track)
    job.collect()
    end_tracking = time.time()
    print(f'{da_particles} tracked in {(end_tracking - start_tracking)/60.:.4f} mins')
    
    parts = job.output.particles[0]
    shape = A1_A2_in_sigma.shape[:2]
    
    init_dict = {'x'    : init_denormalized_6D[:,0].reshape(shape),
                 'px'   : init_denormalized_6D[:,1].reshape(shape),
                 'y'    : init_denormalized_6D[:,2].reshape(shape),
                 'py'   : init_denormalized_6D[:,3].reshape(shape),
                 'tau'  : init_denormalized_6D[:,4].reshape(shape),
                 'ptau' : init_denormalized_6D[:,5].reshape(shape),
                 'A1'   : A1_A2_in_sigma[:,:,0],
                 'A2'   : A1_A2_in_sigma[:,:,1],
                 'J1'   : J1.reshape(shape)/se1**2,
                 'J2'   : J2.reshape(shape)/se2**2
                 }
    
    last_dict = {'x'    : parts.x.reshape(shape),
                 'px'   : parts.px.reshape(shape),
                 'y'   : parts.y.reshape(shape),
                 'py'   : parts.py.reshape(shape),
                 'zeta'   : parts.zeta.reshape(shape),
                 'delta'   : parts.delta.reshape(shape),
                 'at_turn'   : parts.at_turn.reshape(shape)
                }
    
    time_dict = {'setup_time_mins' : (end_setup_time - start_time)/60.,
                 'tracking_time_mins' : (end_tracking - start_tracking)/60.
                }
    
    kfm.dict_to_h5(init_dict, args.output, group='input', readwrite_opts='w')
    kfm.dict_to_h5(last_dict, args.output, group='output', readwrite_opts='a')
    kfm.dict_to_h5(sim_input['partCO'], args.output, group='closed-orbit', readwrite_opts='a')
    kfm.dict_to_h5(optics, args.output, group='optics', readwrite_opts='a')
    kfm.dict_to_h5(time_dict, args.output, group='time', readwrite_opts='a')

    if args.copy_destination is not None:
        shutil.copy(args.output, args.copy_destination + args.output)
elif args.jobtype == 'FMA':
    start_tracking = time.time()
    ecloud_lattice.fma_tracking(distance_between_tunes = 1000, until_turn = 20000, num_stores=n_stores)
    end_tracking = time.time()
    
    output_to_save = {'turn_q'        : np.array(ecloud_lattice.turn_q_list), 
                      'tune_is_valid' : np.array(ecloud_lattice.tune_is_valid_list), 
                      'q1'            : np.array(ecloud_lattice.q1_list),
                      'q2'            : np.array(ecloud_lattice.q2_list),
                      'qx'            : np.array(ecloud_lattice.qx_list),
                      'qy'            : np.array(ecloud_lattice.qy_list),
                      'J1'            : J1/se1**2,
                      'J2'            : J2/se2**2,
                      'setup_time_mins' : (end_setup_time - start_time)/60.,
                      'tracking_time_mins' : (end_tracking - start_tracking)/60.
                     }
    
    kfm.dict_to_h5(output_to_save, args.output)
    if args.copy_destination is not None:
        shutil.copy(args.output, args.copy_destination + args.output)


