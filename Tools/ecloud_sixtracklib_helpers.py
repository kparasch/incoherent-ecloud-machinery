import numpy as np
import time

import sixtracklib
import pysixtrack
import cobjects
import NAFFlib
import kostas_filemanager as kfm

class LatticeWithEclouds:

    def __init__(self, line, eclouds_info, particles_set, device=None):

        self.tricub_data_buffer = cobjects.CBuffer()

        self.tricub_data = {}
        self.tricub_data_indices = {}

        self.tricubs = {}
        self.tricub_indices = {}

        self.tricub_data_buffer_ids = {}

        self.elements = sixtracklib.Elements()
        self.tricub_data_buffer = cobjects.CBuffer()
        self.eclouds_info = eclouds_info

        self.tune_is_valid_list = []
        self.turn_q_list = []
        self.q1_list = []
        self.q2_list = []
        self.qx_list = []
        self.qy_list = []

        self.n_particles = len(particles_set.particles[0].particle_id)
 

        ecloud_list = eclouds_info['length'].keys()

        self.clean_line(line, ecloud_list)
        
        for element, element_name in zip(line.elements, line.element_names):
            element_type = element.__class__.__name__
            if element_name in ecloud_list:
                tc_index = self.elements.cbuffer.n_objects

                tc = sixtracklib.TriCub(cbuffer=self.elements.cbuffer)
                tc.x_shift = self.eclouds_info['x_CO'][element_name]
                tc.y_shift = self.eclouds_info['y_CO'][element_name]
                tc.tau_shift = self.eclouds_info['tau_CO'][element_name]
                tc.length = self.eclouds_info['length'][element_name]

                self.tricubs[element_name] = tc
                self.tricub_indices[element_name] = tc_index
            else:    
                getattr(self.elements, element_type)(**element.to_dict(keepextra=True))

        self.job = sixtracklib.TrackJob(self.elements, particles_set, device=device)

        return

    def clean_line(self, line, ecloud_list):
        
        for ii, elname in enumerate(line.element_names):
            if elname in ecloud_list:
                line.elements[ii] = pysixtrack.elements.Multipole(knl=[100])

        line.remove_inactive_multipoles(inplace=True)
        line.remove_zero_length_drifts(inplace=True)
        line.merge_consecutive_drifts(inplace=True)

        return

    def add_tricub_data(self, fname, data_name, max_z=None):
        grid = kfm.h5_to_obj(fname, group='grid')

        if max_z is not None:
            slices = np.where(abs(grid.zg) < max_z)[0]
            nz = len(slices)
            z0 = grid.zg[slices[0]]
        else:
            nz = grid.Nz
            slices = range(nz)
            z0 = grid.z0
            
        self.tricub_data_indices[data_name] = self.tricub_data_buffer.n_objects
        self.tricub_data[data_name] = sixtracklib.TriCubData(cbuffer=self.tricub_data_buffer,
                                                        nx=grid.Nx, ny=grid.Ny, nz=nz)

        tricub_data = self.tricub_data[data_name]

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

        nx = grid.Nx
        ny = grid.Ny
        for kk in slices:
            kk0 = kk - slices[0]
            if kk0%10 == 0:
                print(f'Loading {data_name} ecloud data ({kk0}/{nz})')
#            start_reading_time = time.time()
            phi = kfm.h5_to_dict(fname, group=f'slices/slice{kk}')['phi'] 
            for ll in range(8):
                phi[:,:,ll] *= scale[ll]
#            end_reading_time = time.time()
            
            index = 8 * (nx * ( ny * kk0))
            tricub_data.table_addr[index: index + 8 * nx * ny] = phi[:, :, :].transpose(1,0,2).flatten()

#            for jj in range(ny):
#                index = 8 * (nx * (jj + ny * kk0))
#                tricub_data.table_addr[index: index + 8 * nx] = phi[:, jj, :].flatten()

##            for ii in range(nx):
##                for jj in range(ny)
##                    for ll in range(8):
##                        tricub_data.table_addr[ll + 8 * (ii + nx * (jj + ny * kk0))] = phi[ii,jj,ll]


#            for ii in range(nx):
#                for jj in range(ny):
#                    for ll in range(8):
#                        if tricub_data.table_addr[ll + 8 * (ii + nx * (jj + ny * kk0))] != phi[ii,jj,ll]:
#                            print('False!')

#            end_writing_time = time.time()
#            print(f'Reading time: {end_reading_time - start_reading_time:.2f}s')
#            print(f'Writing time: {end_writing_time - end_reading_time:.2f}s')
            #breakpoint()
        self.tricub_data_buffer_ids[data_name] = self.job.add_stored_buffer(cbuffer=self.tricub_data_buffer)

        return

    def finalize_assignments(self, tricub_to_tricubdata):
        tricubs = self.tricubs.keys()

        for tc in tricubs:
            tcdata = tricub_to_tricubdata[tc]
            sixtracklib.TriCub_buffer_create_assign_address_item(
                self.job, self.tricub_indices[tc], 
                self.tricub_data_buffer_ids[tcdata], 
                self.tricub_data_indices[tcdata]
            )
        self.job.commit_address_assignments()
        self.job.assign_all_addresses()

        return

    def set_optics_CO(self, optics, partCO):
        self.optics = optics
        self.partCO = pysixtrack.Particles(**partCO)

    def fma_tracking(self, distance_between_tunes=10000, until_turn=20000, num_stores=101):
        stores_center = num_stores//2
        prev_turn_to_track = 0

        last_kk = int(np.ceil(until_turn/distance_between_tunes))

        for kk in range(last_kk+1):
            if kk == 0:
                turn_to_track = num_stores
                if num_stores == distance_between_tunes:
                    continue
            else:
                turn_to_track = np.min([kk*distance_between_tunes, until_turn])
            start_tracking = time.time()
#            print(turn_to_track%num_stores)
            self.job.track_until(turn_to_track)
            end_tracking = time.time()
            time_per_turn = (end_tracking-start_tracking)/(turn_to_track - prev_turn_to_track)
            print(f'Tracked until turn {turn_to_track}/{until_turn}, time/turn = {time_per_turn*1000.:.3f} ms')
            self.job.collect()

            monitor_last_turn = turn_to_track%num_stores
            q1, q2, qx, qy = self.calculate_tunes(self.job.output.particles[0], num_stores, shift=-monitor_last_turn)
            self.turn_q_list.append( turn_to_track - 1 - stores_center)
            mask = self.job.output.particles[0].at_turn.reshape(num_stores, self.n_particles)[-1,:] == turn_to_track - 1
#            print(self.job.output.particles[0].at_turn.reshape(num_stores, self.n_particles)[monitor_last_turn-1,-5:])
#            print(np.roll(self.job.output.particles[0].at_turn.reshape(num_stores, self.n_particles), -monitor_last_turn,axis=0)[-5:,0])
            self.tune_is_valid_list.append(mask)
            self.q1_list.append(q1)
            self.q2_list.append(q2)
            self.qx_list.append(qx)
            self.qy_list.append(qy)
            prev_turn_to_track = turn_to_track


    def calculate_tunes(self, particles, n_turns, shift=0):
        n_particles = self.n_particles
        phys_coords = np.empty([n_turns, n_particles, 6])
        
        phys_coords[:,:,0] = np.roll(particles.x.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.x
        phys_coords[:,:,1] = np.roll(particles.px.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.px
        phys_coords[:,:,2] = np.roll(particles.y.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.y
        phys_coords[:,:,3] = np.roll(particles.py.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.py
        phys_coords[:,:,4] = np.roll(particles.zeta.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.zeta
        phys_coords[:,:,5] = np.roll(particles.delta.reshape(n_turns, n_particles), shift, axis=0) - self.partCO.delta

        norm_coords = np.tensordot( self.optics['invW'], phys_coords, [1,2]).transpose(1,2,0)

        q1 = NAFFlib.multiparticle_tunes(norm_coords[:,:,0].T - 1.j*norm_coords[:,:,1].T).real
        q2 = NAFFlib.multiparticle_tunes(norm_coords[:,:,2].T - 1.j*norm_coords[:,:,3].T).real

        qx = NAFFlib.multiparticle_tunes(phys_coords[:,:,0].T).real
        qy = NAFFlib.multiparticle_tunes(phys_coords[:,:,2].T).real
        
        return q1, q2, qx, qy


