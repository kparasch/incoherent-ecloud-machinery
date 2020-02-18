import pysixtrack
import sixtracklib
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Tools')
import distribution
import ecloud_sixtracklib_helpers as ec_stl
import NAFFlib

from scipy.constants import c as clight

fOptics = 'optics.pkl'
fLine = 'line_with_ecloud_markers_and_collimators.pkl'
fPartCO = 'part_on_CO.pkl'
fEclouds = 'eclouds_info.pkl'

with open(fLine, 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid),keepextra=True)


with open(fOptics, 'rb') as fid:
    optics = pickle.load(fid)

with open(fPartCO, 'rb') as fid:
    partCO = pickle.load(fid)

with open(fEclouds, 'rb') as fid:
    eclouds_info = pickle.load(fid)

ecloud_scale = 1.
for key in eclouds_info['length'].keys():
    eclouds_info['length'][key] *= ecloud_scale/(optics['beta0']*optics['p0c_eV']*clight)

n_turns = 100
n_particles = 1000
n_sigma = 5
ptau_max = 0.#5.e-4#2.7e-4
epsn_1 = 1.7e-6
epsn_2 = 1.7e-6
seed = 10

se1 = np.sqrt(epsn_1/optics['gamma0']/optics['beta0'])
se2 = np.sqrt(epsn_2/optics['gamma0']/optics['beta0'])

line.append_element(pysixtrack.elements.BeamMonitor(num_stores=n_turns),'monitor1')

init_denormalized_6D = distribution.get6D_with_fixed_J3(n_particles=n_particles, 
                                           n_sigma=n_sigma, ptau_max=ptau_max, 
                                           epsn_1=epsn_1, epsn_2=epsn_2, 
                                           optics=optics, seed=seed
                                          )

init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])

ecloud_lattice = ec_stl.LatticeWithEclouds(line, eclouds_info, ps, device=None)

ecloud_lattice.add_tricub_data('/home/kparasch/Storage/MyPinches/refined_Pinch10_MTI1.0_MLI1.0_DTO1.0_DLO1.0.h5', 'drift', max_z=0.15)

tricub_to_tricub_data = {}
for key in eclouds_info['length'].keys():
    tricub_to_tricub_data[key] = 'drift'
ecloud_lattice.finalize_assignments(tricub_to_tricub_data)

ecloud_lattice.job.track_until(n_turns)

ecloud_lattice.job.collect()

x = ecloud_lattice.job.output.particles[0].x.reshape(n_turns, n_particles)
px = ecloud_lattice.job.output.particles[0].px.reshape(n_turns, n_particles)
y = ecloud_lattice.job.output.particles[0].y.reshape(n_turns, n_particles)
py = ecloud_lattice.job.output.particles[0].py.reshape(n_turns, n_particles)
zeta = ecloud_lattice.job.output.particles[0].zeta.reshape(n_turns, n_particles)
delta = ecloud_lattice.job.output.particles[0].delta.reshape(n_turns, n_particles)

denorm_coords = np.empty([n_turns,n_particles,6])
denorm_coords[:,:,0] = x - partCO['x']
denorm_coords[:,:,1] = px - partCO['px']
denorm_coords[:,:,2] = y - partCO['y']
denorm_coords[:,:,3] = py - partCO['py']
denorm_coords[:,:,4] = zeta - partCO['zeta']
denorm_coords[:,:,5] = delta - partCO['delta']
invW = optics['invW']

norm_coords = np.tensordot(invW, denorm_coords, [1,2]).transpose(1,2,0)

q1 = NAFFlib.multiparticle_tunes(denorm_coords[:,:,0].T+1.j*denorm_coords[:,:,1].T).real
q2 = NAFFlib.multiparticle_tunes(denorm_coords[:,:,2].T+1.j*denorm_coords[:,:,3].T).real
print(np.max(zeta))
#plt.figure(1)
#plt.plot(norm_coords[:,:,0]/se1, norm_coords[:,:,1]/se1,'.')
#plt.xlabel('x hat')
#plt.ylabel('px hat')
#plt.figure(2)
#plt.plot(norm_coords[:,:,2]/se2, norm_coords[:,:,3]/se2,'.')
#plt.xlabel('y hat')
#plt.ylabel('py hat')
#plt.figure(3)
#plt.plot(denorm_coords[:,:,4], denorm_coords[:,:,5],'.')
#plt.xlabel('zeta')
#plt.ylabel('delta')
plt.figure(4)
plt.plot(q1,q2,'.')
plt.xlabel('q1')
plt.ylabel('q2')
plt.show()
