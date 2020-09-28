import h5py
import time
import sys
sys.path.append('../Tools')
import kostas_filemanager as kfm
import ecloud_sixtracklib_helpers as ec_stl
import distribution

import pysixtrack
import sixtracklib

import pickle
import argparse
import shutil

parser = argparse.ArgumentParser(description='Long tracking simulations with e-clouds')
parser.add_argument('--ecloud', dest='do_ecloud', action='store_true')
parser.add_argument('--device', nargs='?', default=None, type=str)
parser.add_argument('--copy_destination', nargs='?', default=None, type=str)
parser.add_argument('--ecloud_strength', nargs='?', default=1., type=float)
parser.add_argument('--line_folder', nargs='?', default='./', type=str)
parser.add_argument('--pinch', nargs='?', type=str)
parser.add_argument('--ptau_max', nargs='?', default=0., type=float)
parser.add_argument('--max_tau', nargs='?', default=0.1, type=float)
parser.add_argument('--output', nargs='?', default='temp.h5', type=str)
parser.add_argument('--skip_turns', nargs='?', default=10000, type=int)
parser.add_argument('--turns_per_checkpoint', nargs='?', default=1000000, type=int)
parser.add_argument('--last_checkpoint', nargs='?', default=10, type=int)
parser.add_argument('--n_particles', nargs='?', default=20000, type=int)
parser.add_argument('--n_sigma', nargs='?', default=5.5, type=float)
parser.add_argument('--seed', nargs='?', default=0, type=int)
args = parser.parse_args()
if args.line_folder[-1] != '/':
    args.line_folder += '/'
if args.copy_destination is not None:
    if args.copy_destination[-1] != '/':
        args.copy_destination += '/'

args_dict = vars(args).copy()
for key in args_dict.keys():
    if args_dict[key] is None:
        args_dict[key] = 'None'

line_folder = args.line_folder
device = args.device
ptau_max = args.ptau_max
pinch_path = args.pinch
ecloud_scale = args.ecloud_strength
do_ecloud = args.do_ecloud
max_tau = args.max_tau
output_file = args.output
skip_turns = args.skip_turns
turns_per_checkpoint = args.turns_per_checkpoint
last_checkpoint = args.last_checkpoint
n_sigma = args.n_sigma
copy_destination = args.copy_destination
n_particles = args.n_particles
seed = args.seed

if copy_destination is not None:
    shutil.copy(copy_destination + output_file, output_file)
#pinch_path = 'Pinches/refined_LHC_ArcDip_1.35sey_0.7e11ppb_symm2D_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5'


fOptics = line_folder + 'optics.pkl'
fLine = line_folder + 'line_with_ecloud_markers_and_collimators.pkl'
fPartCO = line_folder + 'part_on_CO.pkl'
fEclouds = line_folder + 'eclouds_info.pkl'
checkpoint = kfm.h5_to_dict(output_file, group='checkpoint')['checkpoint']
print(f'Starting from checkpoint: {checkpoint}!')

with open(fLine, 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid),keepextra=True)

with open(fOptics, 'rb') as fid:
    optics = pickle.load(fid)

with open(fPartCO, 'rb') as fid:
    partCO = pickle.load(fid)

if do_ecloud:
    with open(fEclouds, 'rb') as fid:
        eclouds_info = pickle.load(fid)
    
    for key in eclouds_info['length'].keys():
        eclouds_info['length'][key] *= ecloud_scale/(optics['beta0']*optics['p0c_eV'])
        
    old_optics = optics.copy()
    new_optics = ec_stl.update_optics(line.copy(), eclouds_info.copy(), optics.copy(), partCO.copy(), pinch_path)
    optics = new_optics

start_time = time.time()


n_stores1 = 100
n_stores2 = 1
epsn_1 = 2.0e-6
epsn_2 = 2.0e-6

line.append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores1, start=checkpoint*turns_per_checkpoint + skip_turns-1, skip=skip_turns, is_rolling=True),'monitor1')
line.append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores2, start=checkpoint*turns_per_checkpoint, is_rolling=True),'monitor2')

if checkpoint == 0:
    init_denormalized_6D = distribution.get6D_with_matched_J3_shell(n_particles=n_particles, 
                                                                    n_sigma=n_sigma, 
                                                                    ptau_max=ptau_max, 
                                                                    epsn_1=epsn_1, 
                                                                    epsn_2=epsn_2,
                                                                    optics=optics, 
                                                                    seed=seed
                                                                   )

    init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

    init_dict = {'x': init_denormalized_6D[:,0], 'px': init_denormalized_6D[:,1],
                 'y': init_denormalized_6D[:,2], 'py': init_denormalized_6D[:,3],
                 'tau': init_denormalized_6D[:,4], 'ptau': init_denormalized_6D[:,5]
                }

    ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])
else:
    ps = distribution.load_state(output_file, p0c_eV=optics['p0c_eV'])

if do_ecloud:
    ecloud_lattice = ec_stl.LatticeWithEclouds(line, eclouds_info, ps, device=device)
    ecloud_lattice.set_optics_CO(optics, partCO)

    ecloud_lattice.add_tricub_data(pinch_path, 'dipole', max_z=max_tau)
    ecloud_lattice.remove_dipolar_kicks()

    tricub_to_tricub_data = {}
    for key in eclouds_info['length'].keys():
        tricub_to_tricub_data[key] = 'dipole'
    ecloud_lattice.finalize_assignments(tricub_to_tricub_data)
    job = ecloud_lattice.job
else:
    print(f'Number of elements in line before cleaning: {len(line.elements)}')
    line.remove_inactive_multipoles(inplace=True)
    line.remove_zero_length_drifts(inplace=True)
    line.merge_consecutive_drifts(inplace=True)
    print(f'Number of elements in line after cleaning: {len(line.elements)}')
    elements = sixtracklib.Elements.from_line(line)
    job = sixtracklib.TrackJob(elements, ps, device=device)

end_setup_time = time.time()
print(f'Setting up time: {(end_setup_time - start_time)/60.:.4f}mins')

start_tracking = time.time()
for cc in range(checkpoint+1, last_checkpoint+1):
    checkpoint_dict = {'checkpoint': cc} 
    skip_dicts, last_dict = ec_stl.track_to_checkpoint(job, 
                          n_particles=n_particles, checkpoint=cc, 
                          checkpoint_turns=turns_per_checkpoint, 
                          monitor1_stores=n_stores1, monitor2_stores=n_stores2, 
                          skip_turns=skip_turns
                                                      )


    if cc == 1:
        kfm.dict_to_h5(init_dict, output_file, group='init', readwrite_opts='a')
        kfm.dict_to_h5(optics, output_file, group='optics', readwrite_opts='a')
        kfm.dict_to_h5(partCO, output_file, group='closed-orbit', readwrite_opts='a')
        kfm.dict_to_h5(last_dict, output_file, group='last', readwrite_opts='a')
        time_dict = {'setup_time_mins' : 0., 'tracking_time_mins' : 0. }
        kfm.dict_to_h5(time_dict, output_file, group='time', readwrite_opts='a')
        kfm.dict_to_h5(args_dict, output_file, group='args', readwrite_opts='a')
    else:
        kfm.overwrite(last_dict, output_file, group='last')

    for (turn, skip_dict) in skip_dicts:
        kfm.dict_to_h5(skip_dict, output_file, group=f'turn{turn:d}', readwrite_opts='a')

    kfm.overwrite(checkpoint_dict, output_file, group='checkpoint')

    if copy_destination is not None:
        shutil.copy(output_file, copy_destination + output_file)

    print(f'Reached checkpoint: {cc}/{last_checkpoint}')
end_tracking = time.time()
print(f'Tracking time: {(end_tracking - start_tracking)/60.:.4f}mins')

time_dict = {'setup_time_mins' : (end_setup_time - start_time)/60.,
             'tracking_time_mins' : (end_tracking - start_tracking)/60.
            }

kfm.overwrite(args_dict, output_file, group='args')
kfm.overwrite(time_dict, output_file, group='time')


