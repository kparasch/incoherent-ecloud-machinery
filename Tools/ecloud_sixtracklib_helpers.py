import numpy as np
import time
import copy

import sixtracklib
import pysixtrack
import cobjects
import NAFFlib
import kostas_filemanager as kfm
import normalization

class LatticeWithEclouds:

    def __init__(self, sim_input, particles_set, ecloud_types=None, device=None):
        
        if ecloud_types is None:
            ecloud_types = []

        self.tricub_data = {}
        self.tricub_data_indices = {}

        self.tricubs = {}
        self.tricub_indices = {}

        self.tricub_data_buffer_ids = {}

        self.elements = sixtracklib.Elements()
        self.tricub_data_buffer = cobjects.CBuffer()

        line = sim_input['line']
        self.eclouds_info = sim_input['eclouds_info']

        self.optics = sim_input['optics']
        self.partCO = pysixtrack.Particles(**sim_input['partCO'])
        

        self.tune_is_valid_list = []
        self.turn_q_list = []
        self.q1_list = []
        self.q2_list = []
        self.qx_list = []
        self.qy_list = []

        self.n_particles = len(particles_set.particles[0].particle_id)
 
        self.ecloud_types = ecloud_types

        ecloud_list = [key for key in self.eclouds_info['type'].keys() if self.eclouds_info['type'][key] in ecloud_types]

        print(f'Number of elements in line before cleaning: {len(line.elements)}')
        self.clean_line(line, ecloud_list)
        print(f'Number of elements in line after cleaning: {len(line.elements)}')
        
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

        nx = len(grid.xg)
        ny = len(grid.yg)
            
        self.tricub_data_indices[data_name] = self.tricub_data_buffer.n_objects
        self.tricub_data[data_name] = sixtracklib.TriCubData(cbuffer=self.tricub_data_buffer,
                                                        nx=nx, ny=ny, nz=nz)

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

        for kk in slices:
            kk0 = kk - slices[0]
            if kk0%10 == 0:
                print(f'Loading {data_name} ecloud data ({kk0}/{nz})')
            phi = kfm.h5_to_dict(fname, group=f'slices/slice{kk}')['phi'] 
            for ll in range(8):
                phi[:,:,ll] *= scale[ll]
            
            index = 8 * (nx * ( ny * kk0))
            tricub_data.table_addr[index: index + 8 * nx * ny] = phi[:, :, :].transpose(1,0,2).flatten()

#        self.tricub_data_buffer_ids[data_name] = self.job.add_stored_buffer(cbuffer=self.tricub_data_buffer)

        return

    def remove_dipolar_kicks(self):
        for this_ecloud_type in self.ecloud_types:
            temp_lattice = sixtracklib.Elements()
            sixtracklib.Drift(cbuffer=temp_lattice.cbuffer)
            temp_tc_index = temp_lattice.cbuffer.n_objects
            temp_tc = sixtracklib.TriCub(cbuffer=temp_lattice.cbuffer)
            first_ecloud = [key for key in self.tricubs.keys() if key.split('.')[1] == this_ecloud_type][0]
            print('First ecloud index: {first_ecloud}')
            temp_tc.length = self.tricubs[first_ecloud].length
            temp_tc.x_shift = 0.
            temp_tc.y_shift = 0.
            temp_tc.tau_shift = 0.
            temp_tc.dipolar_kick_px = 0.
            temp_tc.dipolar_kick_py = 0.
            temp_tc.dipolar_kick_ptau = 0.

            temp_ps = particles_set = sixtracklib.ParticlesSet()
            particles  = particles_set.Particles(num_particles=1)

            temp_part = pysixtrack.Particles(p0c=self.partCO.p0c)
            temp_part.x = 0
            temp_part.px = 0
            temp_part.y = 0
            temp_part.py = 0
            temp_part.tau = 0
            temp_part.ptau = 0
            temp_part.state = 1
            temp_part.partid = 0
            temp_part.elemid = 0
            temp_part.turn = 0
            particles.from_pysixtrack(temp_part, 0)

            temp_job = sixtracklib.TrackJob(temp_lattice, temp_ps, device=None)
            temp_tricub_data_buffer_id = temp_job.add_stored_buffer(cbuffer=self.tricub_data_buffer)

            first_tricub_data = list(self.tricub_data.keys())[0]
            sixtracklib.TriCub_buffer_create_assign_address_item(
                temp_job, temp_tc_index, 
                temp_tricub_data_buffer_id, 
                self.tricub_data_indices[first_tricub_data]
            )
            temp_job.commit_address_assignments()
            temp_job.assign_all_addresses()

            temp_job.track_until(1)

            dipolar_kick_px = particles.px[0]
            dipolar_kick_py = particles.py[0]
            dipolar_kick_ptau = particles.ptau[0]
            print(f'{this_ecloud_type} dipolar kicks, px:{dipolar_kick_px}, py:{dipolar_kick_py}, ptau:{dipolar_kick_ptau}')
            #dipolar_kick_px = 0.* particles.px[0]
            #dipolar_kick_py = 0.*particles.py[0]
            #dipolar_kick_ptau = 0.*particles.ptau[0]
            for tc in self.tricubs.keys():
                if tc.split('.')[1] == this_ecloud_type:
                    tc_index = self.tricub_indices[tc]
                    tricub = self.job.beam_elements_buffer.get_object(tc_index)
                    tricub.dipolar_kick_px = dipolar_kick_px
                    tricub.dipolar_kick_py = dipolar_kick_py
                    tricub.dipolar_kick_ptau = dipolar_kick_ptau
            self.job.push_beam_elements()

        return

    def finalize_assignments(self):
        tricub_to_tricubdata = {}
        for key in self.eclouds_info['type'].keys():
            this_type = self.eclouds_info['type'][key]
            if this_type in self.ecloud_types:
                tricub_to_tricubdata[key] = this_type

        tricubs = self.tricubs.keys()

        if tricub_to_tricubdata:
            tricub_data_buffer_id = self.job.add_stored_buffer(cbuffer=self.tricub_data_buffer)
            for tc in tricubs:
                tcdata = tricub_to_tricubdata[tc]
                sixtracklib.TriCub_buffer_create_assign_address_item(
                    self.job, self.tricub_indices[tc], 
                    tricub_data_buffer_id, 
                    self.tricub_data_indices[tcdata]
                )
            self.job.commit_address_assignments()
            self.job.assign_all_addresses()

        return

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


def track_to_checkpoint(job, n_particles, checkpoint=1, checkpoint_turns=1e6, monitor1_stores=100, monitor2_stores=1, skip_turns=10000):
    start_time = time.time()
    turn_to_track = checkpoint * checkpoint_turns
    job.track_until(turn_to_track)
    end_time = time.time()
    print(f'Tracking time (1M): {(end_time - start_time)/60.:.4f}mins')
    start_time = time.time()
    job.collect()


    monitor1 = job.output.particles[0]
    monitor2 = job.output.particles[1]

    p0c = monitor2.p0c[0]
    
    x_skip = monitor1.x.reshape(monitor1_stores, n_particles)
    px_skip = monitor1.px.reshape(monitor1_stores, n_particles)
    y_skip = monitor1.y.reshape(monitor1_stores, n_particles)
    py_skip = monitor1.py.reshape(monitor1_stores, n_particles)
    zeta_skip = monitor1.zeta.reshape(monitor1_stores, n_particles)
    delta_skip = monitor1.delta.reshape(monitor1_stores, n_particles)
    at_turn_skip = monitor1.at_turn.reshape(monitor1_stores, n_particles)

    skip_dicts = []
    for i in range(monitor1_stores):
        zeta = zeta_skip[i]
        delta = delta_skip[i]

        temp_part = pysixtrack.Particles(p0c=p0c)
        temp_part.zeta = zeta
        temp_part.delta = delta

        tau = temp_part.tau
        ptau = temp_part.ptau
        x = x_skip[i]
        px = px_skip[i]
        y = y_skip[i]
        py = py_skip[i]

        supposed_turn = turn_to_track - checkpoint_turns + (i + 1)* skip_turns
        at_turn = at_turn_skip[i] + 1
    #    print(at_turn, supposed_turn)
        mask = at_turn == supposed_turn
        not_mask = np.logical_not(mask)

        x[not_mask] = 0.
        px[not_mask] = 0.
        y[not_mask] = 0.
        py[not_mask] = 0.
        tau[not_mask] = 0.
        ptau[not_mask] = 0.
        at_turn[not_mask] = 0.

        skip_dict = {'x'    : x,
                     'px'   : px,
                     'y'    : y,
                     'py'   : py, 
                     'tau'  : tau,
                     'ptau' : ptau,
                     'at_turn' : at_turn
                    }
        skip_dicts.append((supposed_turn, skip_dict))
    
    x_last = monitor2.x.reshape(monitor2_stores, n_particles)
    px_last = monitor2.px.reshape(monitor2_stores, n_particles)
    y_last = monitor2.y.reshape(monitor2_stores, n_particles)
    py_last = monitor2.py.reshape(monitor2_stores, n_particles)
    zeta_last = monitor2.zeta.reshape(monitor2_stores, n_particles)
    delta_last = monitor2.delta.reshape(monitor2_stores, n_particles)
    at_turn_last = monitor2.at_turn.reshape(monitor2_stores, n_particles)
    state_last = monitor2.state.reshape(monitor2_stores, n_particles)

    temp_part = pysixtrack.Particles(p0c=p0c)
    temp_part.zeta = zeta_last
    temp_part.delta = delta_last
    tau_last = temp_part.tau
    ptau_last = temp_part.ptau

    last_dict = {'x' : x_last.flatten(),
                 'px' : px_last.flatten(),
                 'y'  : y_last.flatten(),
                 'py' : py_last.flatten(),
                 'tau' : tau_last.flatten(),
                 'ptau' : ptau_last.flatten(),
                 'at_turn' : at_turn_last.flatten(),
                 'state' : state_last.flatten()
                }
    

    end_time = time.time()
    print(f'Processing time: {(end_time - start_time)/60.:.4f}mins')
    return skip_dicts, last_dict


def update_optics(sim_input, ecloud_sources, d=1.e-8):
    
    sim_input_copy = copy.deepcopy(sim_input)

    new_optics = sim_input_copy['optics']

    init_part = pysixtrack.Particles(**sim_input_copy['partCO'])
    init_part.x     += np.array([0., 1.*d, 0., 0., 0., 0., 0.])
    init_part.px    += np.array([0., 0., 1.*d, 0., 0., 0., 0.])
    init_part.y     += np.array([0., 0., 0., 1.*d, 0., 0., 0.])
    init_part.py    += np.array([0., 0., 0., 0., 1.*d, 0., 0.])
    init_part.tau   += np.array([0., 0., 0., 0., 0., 1.*d, 0.])
    init_part.ptau  += np.array([0., 0., 0., 0., 0., 0., 1.*d])
    
    n_part = 7

    ps = sixtracklib.ParticlesSet()
    p = ps.Particles(num_particles=n_part)

    for i_part in range(n_part):
        part = pysixtrack.Particles(p0c=init_part.p0c)

        part.x    = init_part.x[i_part]
        part.px   = init_part.px[i_part]
        part.y    = init_part.y[i_part]
        part.py   = init_part.py[i_part]
        part.tau  = init_part.tau[i_part]
        part.ptau = init_part.ptau[i_part]

        part.partid = i_part
        part.state  = 1
        part.elemid = 0
        part.turn   = 0

        p.from_pysixtrack(part, i_part)

    ecloud_lattice = LatticeWithEclouds(sim_input_copy, ps, ecloud_types=list(ecloud_sources.keys()), device=None)
    
    for key in ecloud_sources.keys():
        ecloud_lattice.add_tricub_data(ecloud_sources[key], key, max_z=0.05)
    ecloud_lattice.remove_dipolar_kicks()
    ecloud_lattice.finalize_assignments()
    job = ecloud_lattice.job

    job.track_until(1)

    X_init = np.empty([6,7])
    X_fin = np.empty([6,7])
    X_init[0,:] = init_part.x
    X_init[1,:] = init_part.px
    X_init[2,:] = init_part.y
    X_init[3,:] = init_part.py
    X_init[4,:] = init_part.tau
    X_init[5,:] = init_part.ptau

    fin_part = pysixtrack.Particles(p0c=init_part.p0c)
    fin_part.x     = p.x.flatten()
    fin_part.px    = p.px.flatten()
    fin_part.y     = p.y.flatten()
    fin_part.py    = p.py.flatten()
    fin_part.zeta  = p.zeta.flatten()
    fin_part.delta = p.delta.flatten()

    X_fin[0,:] = fin_part.x
    X_fin[1,:] = fin_part.px
    X_fin[2,:] = fin_part.y
    X_fin[3,:] = fin_part.py
    X_fin[4,:] = fin_part.tau
    X_fin[5,:] = fin_part.ptau

    m = X_fin[:, 0] - X_init[:, 0]
    M = np.empty([6, 6])
    for j in range(6):
        M[:,j] = (X_fin[:,j+1] - X_fin[:,0]) / d

    Ms = normalization.healy_symplectify(M)
    W, invW, R = normalization.linear_normal_form(Ms)
    q1 = np.arccos(R[0,0])/np.pi/2.
    q2 = np.arccos(R[2,2])/np.pi/2.
    Qs = np.arccos(R[4,4])/np.pi/2.
    
    new_optics['old_W'] = copy.deepcopy(new_optics['W'])
    new_optics['old_invW'] = copy.deepcopy(new_optics['invW'])
    new_optics['old_R'] = copy.deepcopy(new_optics['R'])

    new_optics['W'] = W
    new_optics['invW'] = invW
    new_optics['R'] = R
    new_optics['new_Qs'] = Qs
    new_optics['new_q1'] = q1
    new_optics['new_q2'] = q2

    int_q1 = int(new_optics['q1'])
    int_q2 = int(new_optics['q2'])
    old_q1 = new_optics['q1'] - int_q1
    old_q2 = new_optics['q2'] - int_q2
    old_q3 = new_optics['Qs']
    new_q1 = new_optics['new_q1']
    new_q2 = new_optics['new_q2']
    new_q3 = new_optics['new_Qs']
    print(f'Int. tunes: int(q1) = {int_q1}, int(q2) = {int_q2}')
    print(f'Old tunes: q1 = {old_q1:.4f}, q2 = {old_q2:.4f}, q3 = {old_q3:.6f}')
    print(f'New tunes: q1 = {new_q1:.4f}, q2 = {new_q2:.4f}, q3 = {new_q3:.6f}')
    print(f'Tune shift: \u0394q1 = {new_q1-old_q1:.4f}, \u0394q2 = {new_q2-old_q2:.4f}, \u0394q3 = {new_q3-old_q3:.6f}')

    return new_optics



