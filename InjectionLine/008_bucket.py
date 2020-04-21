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


from scipy.constants import c as clight

start_time = time.time()

line_folder = 'Lines/Line_IMO_0/'
device = 'opencl:0.0'
ptau_max = 6.6e-4#5.e-4#2.7e-4
seed = 10
pinch_path = 'Pinches/refined_Pinch32_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5'
ecloud_scale = 0.1
do_ecloud = True
max_z = 0.22
bucket_acceptance = 1.e-3 #(1.e-3 for injection, 0.36e-3 for collisions)
n_turns = 30
n_stores = n_turns

n_particles = 100
n_particles = 20000



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
        eclouds_info['length'][key] *= ecloud_scale/(optics['beta0']*optics['p0c_eV'])

n_sigma = 5
epsn_1 = 3.5e-6
epsn_2 = 3.5e-6

se1 = np.sqrt(epsn_1/optics['gamma0']/optics['beta0'])
se2 = np.sqrt(epsn_2/optics['gamma0']/optics['beta0'])

line.append_element(pysixtrack.elements.BeamMonitor(num_stores=n_stores,is_rolling=True),'monitor1')

init_denormalized_6D = np.zeros([n_particles,6])
init_denormalized_6D[:,5] = np.linspace(0, ptau_max, n_particles) 

init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])

J1, J2 = distribution.J1_J2_from_physical(init_denormalized_6D, optics['invW'], partCO)

if do_ecloud:
    ecloud_lattice = ec_stl.LatticeWithEclouds(line, eclouds_info, ps, device=device)
    ecloud_lattice.set_optics_CO(optics, partCO)

    ecloud_lattice.add_tricub_data(pinch_path, 'drift', max_z=max_z)
#ecloud_lattice.add_tricub_data('refined_Pinch10_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5', 'drift', max_z=0.15)

    tricub_to_tricub_data = {}
    for key in eclouds_info['length'].keys():
        tricub_to_tricub_data[key] = 'drift'
    ecloud_lattice.finalize_assignments(tricub_to_tricub_data)
    job = ecloud_lattice.job
else:
    elements = sixtracklib.Elements.from_line(line)
#    elements.from_line(line)
    job = sixtracklib.TrackJob(elements, ps, device=device)


end_setup_time = time.time()
print(f'Setting up time: {(end_setup_time - start_time)/60.}mins')

start_tracking = time.time()
job.track_until(n_turns)
job.collect()
end_tracking = time.time()

particles = job.output.particles[0]
zeta = particles.zeta.reshape(n_stores,n_particles)
delta = particles.delta.reshape(n_stores,n_particles)
plt.plot(zeta,delta,'.')
plt.xlim(-0.5,0.5)
#ecloud_lattice.fma_tracking(distance_between_tunes = 1000, until_turn = 20000, num_stores=n_stores)

#output_to_save = {'turn_q'        : np.array(ecloud_lattice.turn_q_list), 
#                  'tune_is_valid' : np.array(ecloud_lattice.tune_is_valid_list), 
#                  'q1'            : np.array(ecloud_lattice.q1_list),
#                  'q2'            : np.array(ecloud_lattice.q2_list),
#                  'qx'            : np.array(ecloud_lattice.qx_list),
#                  'qy'            : np.array(ecloud_lattice.qy_list),
#                  'J1'            : J1/se1**2,
#                  'J2'            : J2/se2**2
#                 }

init_dict = {'x'    : init_denormalized_6D[:,0],
             'px'   : init_denormalized_6D[:,1],
             'y'    : init_denormalized_6D[:,2],
             'py'   : init_denormalized_6D[:,3],
             'tau'  : init_denormalized_6D[:,4],
             'ptau' : init_denormalized_6D[:,5],
             'J1'   : J1/se1**2,
             'J2'   : J2/se2**2
             }

#kfm.dict_to_h5(output_to_save, f'fma_tunes_ec_scale{ecloud_scale:.2f}_IMO_0_'+sys.argv[2]+'.h5')
plt.show()
