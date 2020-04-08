import pysixtrack
import sixtracklib
import numpy as np
from scipy.constants import c

def get_fma_distribution(n_particles, n_sigma, ptau_max, epsn_1, epsn_2, optics, seed=0):
    W = optics['W']
    invW = optics['invW']
    gamma0 = optics['gamma0']
    beta0 = optics['beta0']

    e1 = epsn_1/(beta0*gamma0)
    e2 = epsn_2/(beta0*gamma0)

    init_denormalized_6D = np.empty([n_particles, 6])
    init_denormalized_6D[:,0] = 0.
    init_denormalized_6D[:,1] = 0.
    init_denormalized_6D[:,2] = 0.
    init_denormalized_6D[:,3] = 0.
    init_denormalized_6D[:,4] = 0.
    init_denormalized_6D[:,5] = ptau_max
    init_normalized_6D_temp = np.tensordot(invW, init_denormalized_6D, [1,1]).T
    init_normalized_coordinates_sigma = random_full_hypersphere(n_sigma, n_particles, dim=2, seed=seed)

    init_normalized_6D = np.empty([n_particles, 6])
    init_normalized_6D[:,0] = init_normalized_coordinates_sigma[:,0] * np.sqrt(e1)
    init_normalized_6D[:,1] = 0.
    init_normalized_6D[:,2] = init_normalized_coordinates_sigma[:,1] * np.sqrt(e2)
    init_normalized_6D[:,3] = 0.
    init_normalized_6D[:,4] = 0.
    init_normalized_6D[:,5] = init_normalized_6D_temp[:,5]

    init_denormalized_coordinates = np.tensordot(W, init_normalized_6D, [1,1]).T

    return init_denormalized_coordinates


def get6D_with_fixed_J3(n_particles, n_sigma, ptau_max, epsn_1, epsn_2, optics, seed=0):

    alfa = optics['alfa']
    V0 = optics['rf_volt_V']
    freq = optics['rf_freq_Hz']
    gamma0 = optics['gamma0']
    beta0 = optics['beta0']
    p0c_eV = optics['p0c_eV']
    length = optics['length']

    W = optics['W']
    invW = optics['invW']

    eta = alfa - 1./gamma0**2
    A = np.sqrt(2*np.pi*freq*p0c_eV/c*length/4./V0)
    b = 2*np.pi*freq/2./c

    e1 = epsn_1/(beta0*gamma0)
    e2 = epsn_2/(beta0*gamma0)

    phi3 = np.random.uniform(-np.pi, np.pi, n_particles)
    J3 = eta*ptau_max**2
    tau_tilde = np.sqrt(2*J3)*np.cos(phi3)
    ptau_tilde = np.sqrt(2*J3)*np.sin(phi3)
    tau = np.arcsin(tau_tilde*A/np.sqrt(2.))/b
    ptau = ptau_tilde/np.sqrt(eta)/np.sqrt(2.)

    init_denormalized_6D = np.empty([n_particles, 6])
    init_denormalized_6D[:,0] = 0.
    init_denormalized_6D[:,1] = 0.
    init_denormalized_6D[:,2] = 0.
    init_denormalized_6D[:,3] = 0.
    init_denormalized_6D[:,4] = tau
    init_denormalized_6D[:,5] = ptau
    init_normalized_6D_temp = np.tensordot(invW, init_denormalized_6D, [1,1]).T
    init_normalized_coordinates_sigma = random_full_hypersphere(n_sigma, n_particles, dim=4, seed=seed)

    init_normalized_6D = np.empty([n_particles, 6])
    init_normalized_6D[:,0] = init_normalized_coordinates_sigma[:,0] * np.sqrt(e1)
    init_normalized_6D[:,1] = init_normalized_coordinates_sigma[:,1] * np.sqrt(e1)
    init_normalized_6D[:,2] = init_normalized_coordinates_sigma[:,2] * np.sqrt(e2)
    init_normalized_6D[:,3] = init_normalized_coordinates_sigma[:,3] * np.sqrt(e2)
    init_normalized_6D[:,4] = init_normalized_6D_temp[:,4]
    init_normalized_6D[:,5] = init_normalized_6D_temp[:,5]

    init_denormalized_coordinates = np.tensordot(W, init_normalized_6D, [1,1]).T

    return init_denormalized_coordinates

def random_full_hypersphere(sigma, n_particles, dim, seed = 0):
    np.random.seed(seed)
    np.random.uniform((1,10000))
    i = 0
    sphere = np.zeros((n_particles, dim))
    while i < n_particles:
        point = np.random.uniform(low=-sigma, high=sigma, size=(dim,))
        if np.linalg.norm(point) <= sigma:
            sphere[i] = point
            i += 1
    return sphere

def apply_closed_orbit(init_denormalized_coordinates, partCO):
    part = pysixtrack.Particles(**partCO)
    init_denormalized_coordinates[:,0] += part.x
    init_denormalized_coordinates[:,1] += part.px
    init_denormalized_coordinates[:,2] += part.y
    init_denormalized_coordinates[:,3] += part.py
    init_denormalized_coordinates[:,4] += part.tau
    init_denormalized_coordinates[:,5] += part.ptau

    return init_denormalized_coordinates

def J1_J2_from_physical(denormalized_coordinates, invW, partCO):
    part = pysixtrack.Particles(**partCO)
    coords = denormalized_coordinates.copy()
    coords[:,0] -= part.x
    coords[:,1] -= part.px
    coords[:,2] -= part.y
    coords[:,3] -= part.py
    coords[:,4] -= part.zeta
    coords[:,5] -= part.delta
   
    normalized_coords = np.tensordot(invW, coords, [1,1]).T
    J1 = 0.5*(normalized_coords[:,0]**2 + normalized_coords[:,1]**2)
    J2 = 0.5*(normalized_coords[:,2]**2 + normalized_coords[:,3]**2)
    
    return J1, J2
    

def get_sixtracklib_particle_set(init_denormalized_coordinates, p0c_eV):
    n_part = init_denormalized_coordinates.shape[0]

    ps = sixtracklib.ParticlesSet()
    p = ps.Particles(num_particles=n_part)

    for i_part in range(n_part):
        part = pysixtrack.Particles(p0c=p0c_eV)

        part.x    = init_denormalized_coordinates[i_part, 0]
        part.px   = init_denormalized_coordinates[i_part, 1]
        part.y    = init_denormalized_coordinates[i_part, 2]
        part.py   = init_denormalized_coordinates[i_part, 3]
        part.tau  = init_denormalized_coordinates[i_part, 4]
        part.ptau = init_denormalized_coordinates[i_part, 5]

        part.partid = i_part
        part.state  = 1
        part.elemid = 0
        part.turn   = 0
        
        p.from_pysixtrack(part, i_part)

    return ps

