import numpy as np
import sixtracklib
import pysixtrack
import cobjects
import kostas_filemanager as kfm

def push_tricub_data(fname, data_name, max_z=None, device='opencl:0.0', p0c=1., length=1.):
    lattice = sixtracklib.Elements()
    sixtracklib.Drift(cbuffer=lattice.cbuffer)
    tc_index = lattice.cbuffer.n_objects
    tc = sixtracklib.TriCub(cbuffer=lattice.cbuffer)
    tc.length = length 
    tc.x_shift = 0.
    tc.y_shift = 0.
    tc.tau_shift = 0.
    tc.dipolar_kick_px = 0.
    tc.dipolar_kick_py = 0.
    tc.dipolar_kick_ptau = 0.

    particles_set = sixtracklib.ParticlesSet()
    particles = particles_set.Particles(num_particles=1)
    part = pysixtrack.Particles(p0c=p0c)
    part.x = 0.
    part.px = 0.
    part.y = 0.
    part.py = 0.
    part.tau = 0.
    part.ptau = 0.
    part.state = 1
    part.partid = 0
    part.elemid = 0
    part.turn = 0
    particles.from_pysixtrack(part, 0)

    job = sixtracklib.TrackJob(lattice, particles_set, device=device)

    grid  = kfm.h5_to_obj(fname, group='grid')

    if max_z is not None:
        slices = np.where(abs(grid.zg) < max_z)[0]
        nz = len(slices)
        z0 = grid.zg[slices[0]]
    else:
        nz = grid.Nz
        slices = range(nz)
        z0 = grid.z0

    nx = len(grid.xg)
    ny = len(grid.yg)
        
    tricub_data_buffer = cobjects.CBuffer()
    tricub_data_index = tricub_data_buffer.n_objects
    tricub_data = sixtracklib.TriCubData(cbuffer=tricub_data_buffer,
                                          nx=nx, ny=ny, nz=nz)

    tricub_data.x0 = grid.x0
    tricub_data.y0 = grid.y0
    tricub_data.z0 = z0

    tricub_data.dx = grid.dx
    tricub_data.dy = grid.dy
    tricub_data.dz = grid.dz


    if kfm.h5_to_dict(fname, group='settings')['symmetric2D']:
        tricub_data.mirror_x = 1
        tricub_data.mirror_y = 1
        tricub_data.mirror_z = 0
    else:
        tricub_data.mirror_x = 0
        tricub_data.mirror_y = 0
        tricub_data.mirror_z = 0

    scale = [1., grid.dx, grid.dy, grid.dz, grid.dx * grid.dy, grid.dx * grid.dz,
             grid.dy * grid.dz, (grid.dx * grid.dy) * grid.dz]

    for kk in slices:
        kk0 = kk - slices[0]
        if kk0%10 == 0:
            print(f'Loading {data_name} ecloud data ({kk0}/{nz})')
        phi = kfm.h5_to_dict(fname, group=f'slices/slice{kk}')['phi'] 
        for ll in range(8):
            phi[:,:,ll] *= scale[ll]
        
        index = 8 * (nx * ( ny * kk0))
        tricub_data.table_addr[index: index + 8 * nx * ny] = phi[:, :, :].transpose(1,0,2).flatten()

    tricub_data_buffer_id = job.add_stored_buffer(cbuffer=tricub_data_buffer)
    sixtracklib.TriCub_buffer_create_assign_address_item(job, tc_index, 
                                                         tricub_data_buffer_id,
                                                         tricub_data_index
                                                        )
    job.commit_address_assignments()
    job.assign_all_addresses()
    job.track_until(1)
    job.collect()

    dipolar_kick_px = particles.px[0]
    dipolar_kick_py = particles.py[0]
    dipolar_kick_ptau = particles.ptau[0]
    print(dipolar_kick_px, dipolar_kick_py, dipolar_kick_ptau)

    job.collect_beam_elements()
    job.collect_stored_buffer(tricub_data_buffer_id)
    data_address = job.beam_elements_buffer.get_object(tc_index).data_addr

    return job, data_address, [dipolar_kick_px, dipolar_kick_py, dipolar_kick_ptau]

def remove_dipolar_kicks(ecloud_lattice, dipolar_kicks):
    for tc in ecloud_lattice.tricubs.keys():
        tc_index = ecloud_lattice.tricub_indices[tc]
        tricub = ecloud_lattice.job.beam_elements_buffer.get_object(tc_index)
        tc_length = tricub.length
        tc_length = 1.
        tricub.dipolar_kick_px = dipolar_kicks[0]*tc_length
        tricub.dipolar_kick_py = dipolar_kicks[1]*tc_length
        tricub.dipolar_kick_ptau = dipolar_kicks[2]*tc_length

    return 

def set_address(ecloud_lattice, data_address):
    for tc in ecloud_lattice.tricubs.keys():
        tc_index = ecloud_lattice.tricub_indices[tc]
        tricub = ecloud_lattice.job.beam_elements_buffer.get_object(tc_index)
        tricub.data_addr = data_address

    return 


