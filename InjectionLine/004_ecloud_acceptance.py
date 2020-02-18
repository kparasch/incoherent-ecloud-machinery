import pysixtrack
import sixtracklib
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Tools')
import distribution

fOptics = 'optics.pkl'
fLine = 'line_with_ecloud_markers_and_collimators.pkl'
fPartCO = 'part_on_CO.pkl'

with open(fLine, 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid),keepextra=True)

monitors = 0
for i,elname in enumerate(line.element_names):
    if elname[0:6] == 'ecloud':
        line.elements[i] = pysixtrack.elements.BeamMonitor(num_stores=100, is_rolling=True)
        monitors += 1

with open(fOptics, 'rb') as fid:
    optics = pickle.load(fid)

with open(fPartCO, 'rb') as fid:
    partCO = pickle.load(fid)

n_turns = 100
n_particles = 10000
n_sigma = 10
ptau_max = 2.7e-4
epsn_1 = 1.7e-6
epsn_2 = 1.7e-6
seed = 0

init_denormalized_6D = distribution.get6D_with_fixed_J3(n_particles=n_particles, 
                                           n_sigma=n_sigma, ptau_max=ptau_max, 
                                           epsn_1=epsn_1, epsn_2=epsn_2, 
                                           optics=optics, seed=seed
                                          )

init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])


elements = sixtracklib.Elements()
elements.BeamMonitor(num_stores=n_turns)
elements.append_line(line)

xmax = 0
ymax = 0
tmax = 0
job = sixtracklib.TrackJob( elements, ps, device=None)
first_turns = 10
job.track_until(n_turns)
job.collect()
res = job.output
prr
for j in range(1, monitors+1):
    print(j)
    parts = res.particles[j]
    x = parts.x.reshape(n_turns,n_particles)
    y = parts.x.reshape(n_turns,n_particles)
    zeta = parts.zeta.reshape(n_turns,n_particles)
    at_turn = parts.at_turn.reshape(n_turns,n_particles)
    for i in range(first_turns, n_turns):
        mask = at_turn[i,:] == i
        xmax_t = np.max(np.abs(x[i,:][mask]))
        ymax_t = np.max(np.abs(y[i,:][mask]))
        tmax_t = np.max(np.abs(zeta[i,:][mask]))
        if xmax < xmax_t: xmax = xmax_t
        if ymax < ymax_t: ymax = ymax_t
        if tmax < tmax_t: tmax = tmax_t
        print(f'Max( x: {xmax:.6f}, y: {ymax:.6f}, zeta: {tmax:.6f}), {np.sum(mask)} particles survived, Turn: {i}, Mon: {j}')

#x = res.particles[0].x.reshape(n_turns,n_particles)
#px = res.particles[0].px.reshape(n_turns,n_particles)
#y = res.particles[0].y.reshape(n_turns,n_particles)
#py = res.particles[0].py.reshape(n_turns,n_particles)
#zeta = res.particles[0].zeta.reshape(n_turns,n_particles)
#delta = res.particles[0].delta.reshape(n_turns,n_particles)
#state = res.particles[0].state.reshape(n_turns,n_particles)
