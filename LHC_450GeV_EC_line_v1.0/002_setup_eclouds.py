import numpy as np
import matplotlib.pyplot as plt
import pickle
import pysixtrack


import sys
sys.path.append('../Tools/')
import normalization
import kostas_filemanager as kfm

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--noblock', dest='plot_block', action='store_false')
args = parser.parse_args()

inspect_closed_orbit = True

with open('line_with_ecloud_markers.pkl', 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid),keepextra=True)

with open('ecloud_lengths.pkl', 'rb') as fid:
    ecloud_lengths = pickle.load(fid)

with open('optics.pkl', 'rb') as fid:
    optics = pickle.load(fid)

p0c_eV = optics['p0c_eV']
part = pysixtrack.Particles(p0c = p0c_eV)
part_on_CO, M = normalization.get_CO_and_linear_map(part, line, d=1e-7, tol=1e-9)
part_on_CO.partid = 0
part_on_CO.turn = 0
part_on_CO.state = 0

print(f'before symplectifying det(M) = {np.linalg.det(M)}')
Ms = normalization.healy_symplectify(M)
print(f'after symplectifying det(M) = {np.linalg.det(Ms)}')

W, invW, R = normalization.linear_normal_form(Ms)
Qs = np.arccos(R[4,4])/np.pi/2.

optics['W'] = W
optics['invW'] = invW
optics['R'] = R
optics['Qs'] = Qs

with open('optics.pkl', 'wb') as fid:
    pickle.dump(optics, fid)

with open('part_on_CO.pkl', 'wb') as fid:
    pickle.dump(part_on_CO.to_dict(), fid)

with open('ecloud_beta.pkl', 'rb') as fid:
    ecloud_beta = pickle.load(fid)

ecloud_x_CO = {}
ecloud_y_CO = {}
ecloud_tau_CO = {}
ecloud_type = {}

closed_orbit = line.track_elem_by_elem(part_on_CO)

for i, elname in enumerate(line.element_names):
    if elname[0:6] == 'ecloud':
        ecloud_x_CO[elname] = closed_orbit[i].x
        ecloud_y_CO[elname] = closed_orbit[i].y
        ecloud_tau_CO[elname] = closed_orbit[i].tau
        ecloud_type[elname] = elname.split('.')[1]

eclouds_info = {'length' : ecloud_lengths,
               'x_CO' : ecloud_x_CO,
               'y_CO' : ecloud_y_CO,
               'tau_CO' : ecloud_tau_CO,
               'type' : ecloud_type,
               'beta_x' : {key: ecloud_beta[key]['betx'] for key in ecloud_beta.keys()},
               'beta_y' : {key: ecloud_beta[key]['bety'] for key in ecloud_beta.keys()}
               }

with open('eclouds_info.pkl','wb') as fid:
    pickle.dump(eclouds_info, fid)


if inspect_closed_orbit:
    x = []
    px = []
    y = []
    py = []
    tau = []
    ptau = []
    nturns = 100
    for ii in range(nturns):
        if ii%10 == 0 :
            print(f'{ii/nturns}')
        line.track(part_on_CO)
        x.append(part_on_CO.x)
        px.append(part_on_CO.px)
        y.append(part_on_CO.y)
        py.append(part_on_CO.py)
        tau.append(part_on_CO.tau)
        ptau.append(part_on_CO.ptau)
    
    plt.close(1)
    fig = plt.figure(1)
    ax1 =  fig.add_subplot(321)
    ax1.plot(x)
    ax1.set_ylabel('x')
    ax2 =  fig.add_subplot(322)
    ax2.plot(px)
    ax2.set_ylabel('px')
    ax3 =  fig.add_subplot(323)
    ax3.plot(y)
    ax3.set_ylabel('y')
    ax4 =  fig.add_subplot(324)
    ax4.plot(py)
    ax4.set_ylabel('py')
    ax5 =  fig.add_subplot(325)
    ax5.plot(tau)
    ax5.set_ylabel('tau')
    ax6 =  fig.add_subplot(326)
    ax6.set_ylabel('ptau')
    ax6.plot(ptau)
    
    plt.show(block=args.plot_block)
