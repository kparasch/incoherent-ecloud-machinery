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

with open(fOptics, 'rb') as fid:
    optics = pickle.load(fid)

with open(fPartCO, 'rb') as fid:
    partCO = pickle.load(fid)

n_turns = 400
n_particles = 40000
n_sigma = 10
ptau_max = 2.7e-4
epsn_1 = 1.7e-6
epsn_2 = 1.7e-6
seed = 0
device = 'opencl:0.3'
gamma0 = optics['gamma0']
beta0 = optics['beta0']
e1 = epsn_1/(beta0*gamma0)
e2 = epsn_2/(beta0*gamma0)
ref_epsn = 3.5e-6
coll_r = 5.7*np.sqrt(ref_epsn/epsn_1)
disable_sextupoles = True

for ii in range(nn):
    if line.elements[ii].__class__ is pysixtrack.elements.Multipole:
        order = line.elements[ii].order
        if order == 2:
            line.elements[ii].knl[2] = 0.

#init_denormalized_6D = distribution.get6D_with_fixed_J3(n_particles=n_particles, 
init_denormalized_6D = distribution.get_fma_distribution(n_particles=n_particles, 
                                           n_sigma=n_sigma, ptau_max=ptau_max, 
                                           epsn_1=epsn_1, epsn_2=epsn_2, 
                                           optics=optics, seed=seed
                                          )

init_denormalized_6D = distribution.apply_closed_orbit(init_denormalized_6D, partCO)

ps = distribution.get_sixtracklib_particle_set(init_denormalized_6D, p0c_eV=optics['p0c_eV'])


elements = sixtracklib.Elements()
elements.BeamMonitor(num_stores=n_turns,is_rolling=True)
elements.append_line(line)

xmax = 0
ymax = 0
tmax = 0
job = sixtracklib.TrackJob( elements, ps, device=device)
job.track_until(n_turns)
job.collect()
res = job.output

parts = res.particles[0]
x = parts.x.reshape(n_turns,n_particles)
px = parts.px.reshape(n_turns,n_particles)
y = parts.y.reshape(n_turns,n_particles)
py = parts.py.reshape(n_turns,n_particles)
tau = parts.zeta.reshape(n_turns,n_particles)
ptau = parts.delta.reshape(n_turns,n_particles)
at_turn = parts.at_turn.reshape(n_turns,n_particles)

xmax = np.max(np.abs(x), axis=0)/(np.sqrt(optics['betx']*e1))
ymax = np.max(np.abs(y), axis=0)/(np.sqrt(optics['bety']*e2))

#phys_coords = np.empty([n_turns, n_particles,6])
#phys_coords[:,:,0] = x
#phys_coords[:,:,1] = px
#phys_coords[:,:,2] = y
#phys_coords[:,:,3] = py
#phys_coords[:,:,4] = tau
#phys_coords[:,:,5] = ptau
#
#phys_coords = phys_coords.reshape(n_turns*n_particles,6)
#norm_coords = np.tensordot(optics['invW'], phys_coords, [1,1]).T
#norm_coords = norm_coords.reshape(n_turns, n_particles, 6)
#fJ1 = 0.5*(np.power(norm_coords[:,:,0],2) + np.power(norm_coords[:,:,1],2))/e1
#fJ2 = 0.5*(np.power(norm_coords[:,:,2],2) + np.power(norm_coords[:,:,3],2))/e2
#afJ1 = fJ1[-1,:]
#afJ2 = fJ2[-1,:]
#afA1 = np.sqrt(2*afJ1)
#afA2 = np.sqrt(2*afJ2)
init_norm_coords = np.tensordot(optics['invW'], init_denormalized_6D, [1,1]).T

#for ii in range(6):
#    print(np.max(norm_coords[:,:,ii]))
J1 = 0.5*(np.power(init_norm_coords[:,0],2) + np.power(init_norm_coords[:,1],2))/e1
J2 = 0.5*(np.power(init_norm_coords[:,2],2) + np.power(init_norm_coords[:,3],2))/e2
A1 = np.sqrt(2*J1)
A2 = np.sqrt(2*J2)

turn = n_turns-1
particle_turns = at_turn[-1,:]
mask = particle_turns == turn
not_mask = np.logical_not(mask)
print(sum(mask))
plt.close('all')
#plt.plot(A1[not_mask],A2[not_mask],'r.')
#plt.plot(A1[mask],A2[mask],'b.')
plt.plot(xmax[not_mask],ymax[not_mask],'r.')
plt.plot(xmax[mask],ymax[mask],'b.')
#plt.plot(afA1[mask],afA2[mask],'g.')
plt.axhline(coll_r, color='k',linewidth=3.0)
plt.axvline(coll_r, color='k',linewidth=3.0)
x_circ = np.linspace(0,coll_r,1000)
y_circ = coll_r*np.sqrt(1-(x_circ/coll_r)**2)
plt.plot(x_circ,y_circ, color='k', linewidth=3.0)
x_line = np.linspace(-0.5,10,1000)
coll_deg = 126.91
#plt.plot(x_line, np.tan(coll_deg*np.pi/180.)*x_line + coll_r/np.cos(coll_deg*np.pi/180.))
plt.plot(x_line, np.tan(coll_deg*np.pi/180.)*x_line - coll_r/np.cos(coll_deg*np.pi/180.), color='k', linewidth=3.0)
plt.plot(x_line, np.tan((90.-coll_deg)*np.pi/180.)*x_line + coll_r/np.cos((90-coll_deg)*np.pi/180.), color='k', linewidth=2.0)
#plt.plot(x_line, np.tan((90.-coll_deg)*np.pi/180.)*x_line - coll_r/np.cos((180-coll_deg)*np.pi/180.))
plt.xlim(-0.5,10)
plt.ylim(-0.5,10)
plt.xlabel('A1')
plt.ylabel('A2')
plt.show()
#x = res.particles[0].x.reshape(n_turns,n_particles)
#px = res.particles[0].px.reshape(n_turns,n_particles)
#y = res.particles[0].y.reshape(n_turns,n_particles)
#py = res.particles[0].py.reshape(n_turns,n_particles)
#zeta = res.particles[0].zeta.reshape(n_turns,n_particles)
#delta = res.particles[0].delta.reshape(n_turns,n_particles)
#state = res.particles[0].state.reshape(n_turns,n_particles)
