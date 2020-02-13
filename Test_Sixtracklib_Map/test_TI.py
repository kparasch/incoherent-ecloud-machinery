import sixtracklib as st
import pysixtrack
import cobjects
import numpy as np
from TricubicInterpolation import cTricubic as cTI


device = 'opencl:0.0'
#device = None
n_part = 1000

seed = np.random.randint(10000)
#seed = 7489
print('Seed is %d:'%seed)
np.random.seed(seed)

lattice = st.Elements()
tc_index = lattice.cbuffer.n_objects
tc = st.TriCub(cbuffer=lattice.cbuffer)
tc.length = 1.0


particles_set = st.ParticlesSet()
particles = particles_set.Particles(num_particles=n_part)


nx = 5
ny = 7
nz = 9
A = np.random.rand(nx, ny, nz, 8) * 1.0e-3
dx = 0.001
dy = 0.002
dz = 0.003
x0 = -(nx // 2) * dx
y0 = -(ny // 2) * dy
z0 = -(nz // 2) * dz


TI = cTI.Tricubic_Interpolation(
    A=A, dx=dx, dy=dy, dz=dz, x0=x0, y0=y0, z0=z0, method="Exact"
)

test_x = x0 + TI.ix_bound_up * dx * np.random.rand(n_part)
test_y = y0 + TI.iy_bound_up * dy * np.random.rand(n_part)
test_z = z0 + TI.iz_bound_up * dz * np.random.rand(n_part)

# for x,y,z in zip(test_x,test_y,test_z):
#    print(x,y,z)
#    TI.val(x,y,z)

for i_part in range(n_part):
    part = pysixtrack.Particles()
    part.x = test_x[i_part]
    part.y = test_y[i_part]
    part.tau = test_z[i_part]

    part.partid = i_part
    part.state = 1
    part.elemid = 0
    part.turn = 0
    particles.from_pysixtrack(part, i_part)

job = st.TrackJob(lattice, particles_set, device=device)

tricub_data_buffer = cobjects.CBuffer()
tc_data_index = tricub_data_buffer.n_objects
tc_data = st.TriCubData(cbuffer=tricub_data_buffer, nx=nx, ny=ny, nz=nz)

tc_data.x0 = x0
tc_data.y0 = y0
tc_data.z0 = z0
tc_data.dx = dx
tc_data.dy = dy
tc_data.dz = dz
tc_data.mirror_x = 0
tc_data.mirror_y = 0
tc_data.mirror_z = 0
scale = [1., dx, dy, dz, dx*dy, dx*dz, dy*dz, (dx*dy)*dz]
for ii in range(nx):
    for jj in range(ny):
        for kk in range(nz):
            for ll in range(8):
                tc_data.table_addr[ll + 8 * (ii + nx * (jj + ny * kk ))] = A[
                    ii, jj, kk, ll
                ]*scale[ll]

tricub_data_buffer_id = job.add_stored_buffer(cbuffer=tricub_data_buffer)

st.TriCub_buffer_create_assign_address_item(
    job, tc_index, tricub_data_buffer_id, tc_data_index
)

job.commit_address_assignments()
job.assign_all_addresses()
job.track_until(1)
job.collect()

flag = True
for i_part in range(n_part):
    flag = flag and abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[0] - particles.px[i_part]) < 1.e-12
    flag = flag and abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[1] - particles.py[i_part]) < 1.e-12
    if np.isnan(particles.ptau[i_part]):
        print('ptau at i_part = %d is nan'%i_part)
        if TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[2] < (1 - part.gamma0)/(part.gamma0*part.beta0):
            print('... as it should be.')
        else:
            print('... and should not be.')
            flag = False
    else:
        flag = flag and abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[2] - particles.ptau[i_part]) < 1.e-12
    if not flag:
        print('i_part: ',i_part)
        print(abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[0] - particles.px[i_part]))
        print(abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[1] - particles.py[i_part]))
        print(abs(TI.kick(test_x[i_part], test_y[i_part], test_z[i_part])[2] - particles.ptau[i_part]))
        break

if flag:
    print('Test passed')
else:
    print('Test failed')
#print(flag)
#print(TI.kick(test_x[0], test_y[0], test_z[0]))
#print(particles)
