import pysixtrack
import sixtracklib
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Tools')
import distribution
import ecloud_sixtracklib_helpers as ec_stl
import kostas_filemanager as kfm

import NAFFlib
import time

import argparse
from scipy.constants import c as clight

parser = argparse.ArgumentParser(description='DA simulations with e-clouds')
parser.add_argument('--ecloud', dest='do_ecloud', action='store_true')
#parser.add_argument('--ecloud', dest='do_ecloud', nargs='?', const=True, default=False)
parser.add_argument('--device', nargs='?', default=None, type=str)
parser.add_argument('--ecloud_strength', nargs='?', default=1., type=float)
parser.add_argument('--line_folder', nargs='?', type=str)
parser.add_argument('--pinch', nargs='?', type=str)
parser.add_argument('--ptau_max', nargs='?', default=0., type=float)
parser.add_argument('--max_tau', nargs='?', default=0.1, type=float)
parser.add_argument('--output', nargs='?', default='temp.h5', type=str)
args = parser.parse_args()
if args.line_folder[-1] != '/':
    args.line_folder += '/'

start_time = time.time()

line_folder = args.line_folder
device = args.device
ptau_max = args.ptau_max#5.e-4#2.7e-4
pinch_path = args.pinch
ecloud_scale = args.ecloud_strength
do_ecloud = args.do_ecloud
max_tau = args.max_tau
output_file = args.output



fOptics = line_folder + 'optics.pkl'
fLine = line_folder + 'line_with_ecloud_markers_and_collimators.pkl'
fPartCO = line_folder + 'part_on_CO.pkl'
fEclouds = line_folder + 'eclouds_info.pkl'

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
        eclouds_info['length'][key] *= ecloud_scale/(optics['beta0']*optics['p0c_eV'])/3000.

n_stores = 1
turn_to_track=1000000
n_particles_approx = 20000
n_sigma = 5.7
epsn_1 = 3.5e-6
epsn_2 = 3.5e-6

se1 = np.sqrt(epsn_1/optics['gamma0']/optics['beta0'])
se2 = np.sqrt(epsn_2/optics['gamma0']/optics['beta0'])

line.append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores,is_rolling=True),'monitor1')

init_denormalized_6D, A1_A2_in_sigma, n_particles = distribution.get_DA_distribution(
                                                     n_particles_approx=n_particles_approx, 
                                                     n_sigma=n_sigma, ptau_max=ptau_max, 
                                                     epsn_1=epsn_1, epsn_2=epsn_2, 
                                                     optics=optics
                                                                                     )

init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])

J1, J2 = distribution.J1_J2_from_physical(init_denormalized_6D, optics['invW'], partCO)

if do_ecloud:
    ecloud_lattice = ec_stl.LatticeWithEclouds(line, eclouds_info, ps, device=device)
    ecloud_lattice.set_optics_CO(optics, partCO)

    ecloud_lattice.add_tricub_data(pinch_path, 'drift', max_z=max_tau)

    tricub_to_tricub_data = {}
    for key in eclouds_info['length'].keys():
        tricub_to_tricub_data[key] = 'drift'
    ecloud_lattice.finalize_assignments(tricub_to_tricub_data)
    job = ecloud_lattice.job
else:
    line.remove_inactive_multipoles(inplace=True)
    line.remove_zero_length_drifts(inplace=True)
    line.merge_consecutive_drifts(inplace=True)
    elements = sixtracklib.Elements.from_line(line)
    job = sixtracklib.TrackJob(elements, ps, device=device)


end_setup_time = time.time()
print(f'Setting up time: {(end_setup_time - start_time)/60.}mins')

start_tracking = time.time()
job.track_until(turn_to_track)
job.collect()
end_tracking = time.time()
print(f'Tracking time: {(end_tracking - start_tracking)/60.}mins')

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

kfm.dict_to_h5(init_dict, output_file, group='input', readwrite_opts='w')
kfm.dict_to_h5(last_dict, output_file, group='output', readwrite_opts='a')
kfm.dict_to_h5(partCO, output_file, group='closed-orbit', readwrite_opts='a')
kfm.dict_to_h5(optics, output_file, group='optics', readwrite_opts='a')

args_dict = vars(args)
for key in args_dict.keys():
    if args_dict[key] is None:
        args_dict[key] = 'None'
kfm.dict_to_h5(args_dict, output_file, group='args', readwrite_opts='a')
kfm.dict_to_h5(time_dict, output_file, group='time', readwrite_opts='a')
