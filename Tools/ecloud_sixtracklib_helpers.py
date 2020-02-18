import numpy as np
import time

import sixtracklib
import pysixtrack
import cobjects
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

        self.job = sixtracklib.TrackJob(self.elements, particles_set)

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


