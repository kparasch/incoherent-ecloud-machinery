import numpy as np
import matplotlib.pyplot as plt
import h5py

def plot_init(coll_r, axis):
    x_circ = np.linspace(0,coll_r,1000)
    y_circ = np.sqrt(1 - (x_circ/coll_r)**2)*coll_r
    axis.plot(x_circ, y_circ, 'k:', label='Initial cond.')

def plot_collimation(coll_r, axis, label=None):
    coll_deg = 126.91
    coll_deg = 135
    a = np.tan(coll_deg*np.pi/180.)
    cs = np.abs(np.cos(coll_deg*np.pi/180.))

    x1 = coll_r*(1.-1./cs)/a
    x2 = coll_r
    y1 = a*x1+coll_r/cs
    y2 = a*x2+coll_r/cs


    axis.plot([0,x1],[coll_r,coll_r], color='k', linewidth=3.0)
    axis.plot([x1,x2],[y1,y2], color='k', linewidth=3.0)
    axis.plot([x2,x2],[0,y2], color='k', linewidth=3.0, label='Phys. Aperture')

def plot_DA(fname, axis, label=None, fmt=None, linewidth=2.0):
    #fname = 'DA_IMO_40_ptau_6.6e-4_ecloud.h5'
    
    h5file = h5py.File(fname,'r')
    
    A1 = h5file['input/A1'][()]
    A2 = h5file['input/A2'][()]
    last_turn = h5file['output/at_turn'][()] + 1
    
    mask = last_turn == 1000000
    not_mask = np.logical_not(mask)
    da_mask = np.argmax(not_mask,axis=1)
    
    nangles = A1.shape[0]
    A1_DA = np.empty(nangles)
    A2_DA = np.empty(nangles)
    for ii in range(nangles):
        A1_DA[ii] = A1[ii, da_mask[ii]-1]
        A2_DA[ii] = A2[ii, da_mask[ii]-1]
    
    axis.plot(A1_DA, A2_DA, fmt, label=label)



def plot_DA_blue_red(fname, axis):
    #fname = 'DA_IMO_40_ptau_6.6e-4_ecloud.h5'
    
    h5file = h5py.File(fname,'r')
    
    A1 = h5file['input/A1'][()]
    A2 = h5file['input/A2'][()]
    last_turn = h5file['output/at_turn'][()] + 1
    
    mask = last_turn == 1000000
    not_mask = np.logical_not(mask)
    
    axis.plot(A1[mask], A2[mask], 'b.')
    axis.plot(A1[not_mask], A2[not_mask], 'r.')
