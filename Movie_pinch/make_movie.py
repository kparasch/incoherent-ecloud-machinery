import numpy as np
import matplotlib
import shutil

import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'jet'
import os
plt.style.use('kostas')

from scipy.constants import e as qe

import sys
sys.path.append('../Tools')

import kostas_filemanager as kfm

pinch_in = 'Pinch29.h5'
tag = pinch_in.split('.')[0]
N_discard = 10

folder_movie = 'temp_frames'
os.mkdir(folder_movie)

grid = kfm.h5_to_obj(pinch_in,group='grid')
shapes = [len(grid.xg), len(grid.yg), len(grid.zg)]
ob0 = kfm.h5_to_obj(pinch_in, group='slices/slice0')

phi=np.zeros(shapes)
rho=np.zeros(shapes)
ex=np.zeros(shapes)
ey=np.zeros(shapes)
ez=np.zeros(shapes)

dx = grid.xg[1]-grid.xg[0]
dy = grid.yg[1]-grid.yg[0]
dz = grid.zg[1]-grid.zg[0]

for ii in range(shapes[2]):
#for ii in range(10):
    if ii%20 == 0:
        print(f'Reading {ii}/{shapes[2]}')
    ob = kfm.h5_to_obj(pinch_in, group=f'slices/slice{ii}')
    phi[:,:,ii] = ob.phi[:,:]
    rho[:,:,ii] = ob.rho[:,:]/(-qe)

ex[1:-1,:,:] = -(phi[2:,:,:] - phi[:-2,:,:])/dx
ey[:,1:-1,:] = -(phi[:,2:,:] - phi[:,:-2,:])/dy
ez[:,:,1:-1] = -(phi[:,:,2:] - phi[:,:,:-2])/dz

if N_discard > 0:
    grid.xg = grid.xg[N_discard:-N_discard]
    grid.yg = grid.yg[N_discard:-N_discard]
    grid.zg = grid.zg[N_discard:-N_discard]
    phi = phi[N_discard:-N_discard, N_discard:-N_discard, N_discard:-N_discard]
    rho = rho[N_discard:-N_discard, N_discard:-N_discard, N_discard:-N_discard]
    ex = ex[N_discard:-N_discard, N_discard:-N_discard, N_discard:-N_discard]
    ey = ey[N_discard:-N_discard, N_discard:-N_discard, N_discard:-N_discard]
    ez = ez[N_discard:-N_discard, N_discard:-N_discard, N_discard:-N_discard]

x_obs = 0
y_obs = 0
ix_obs = np.argmin(np.abs(grid.xg-x_obs))
iy_obs = np.argmin(np.abs(grid.yg-y_obs))
plt.close('all')

fig1 = plt.figure(1, figsize=(8, 6*1.5))

ax1 = fig1.add_subplot(3,1,1)
ax2 = fig1.add_subplot(3,1,2, sharex=ax1)
ax3 = fig1.add_subplot(3,1,3, sharex=ax1)

t = np.linspace(0, 2*np.pi,100)

mbl = ax1.pcolormesh(grid.zg,grid.yg,rho[ix_obs, :, :])
plt.colorbar(mappable=mbl, ax=ax1, aspect=5)
mbl = ax2.pcolormesh(grid.zg,grid.yg,ey[ix_obs, :, :])
plt.colorbar(mappable=mbl, ax=ax2, aspect=5)
mbl = ax3.pcolormesh(grid.zg,grid.yg,phi[ix_obs, :, :])
plt.colorbar(mappable=mbl, ax=ax3, aspect=5)

for ax in ax1,ax2,ax3:
    for nn in 1,2,3,4:
        ax.plot(nn*grid.sigma_z_beam*np.cos(t), nn*grid.sigma_y_beam*np.sin(t), color='lightgray', linestyle='-')

z_movie = grid.zg[::-1]

fig2 = plt.figure(2,figsize=(8*1.5, 6*1.5))
fig3 = plt.figure(3,figsize=(9, 5))

rho_max = np.max(rho)
exy_max = np.max([np.max(np.abs(ex)), np.max(np.abs(ey))])
xy_max = np.max([np.max(grid.xg), np.max(grid.yg)])
ez_max = np.max(ez)
ez_min = np.min(ez)
phi_min = np.min(phi)
phi_max = np.max(phi)

rho_min = np.power(np.log10(rho_max)-3, 10)
rho[rho < rho_min] = rho_min

for i_frame, z_obs in enumerate(z_movie):
    print(f'Frame {i_frame}/{len(z_movie)}')

    fig2.clear()
    iz_obs = len(z_movie)-i_frame-1

    axc1 = plt.subplot2grid(shape=(3, 3), loc=(0,0), colspan=2, fig=fig2)
    axc2 = plt.subplot2grid(shape=(3, 3), loc=(1,0), colspan=2, fig=fig2, sharex=axc1)
    axc3 = plt.subplot2grid(shape=(3, 3), loc=(2,0), colspan=2, fig=fig2, sharex=axc1)

    axc1.plot(grid.yg*1e3, rho[ix_obs, :, iz_obs], '.-b', label='$\\mathbf{x=0}$')
    axc1.plot(grid.xg*1e3, rho[:, iy_obs, iz_obs], '.-r', label='$\\mathbf{y=0}$')
    axc1.set_ylim(0, rho_max)
    axc1.set_xlim(-xy_max, xy_max)
    axc1.set_ylabel('$\\mathbf{\\rho\ \ [e^-/m^3]}$')
    axc1.legend(loc='upper right', fontsize=20)

    axc2.plot(grid.yg*1e3, ey[ix_obs, :, iz_obs], '.-b')
    axc2.plot(grid.xg*1e3, ex[:, iy_obs, iz_obs], '.-r')
    axc2.set_ylim(-exy_max, exy_max)
    axc2.set_xlim(-xy_max, xy_max)
    axc2.set_ylabel('$\\mathbf{E_{x,y}\ \ [V/m]}$')

    axc3.plot(grid.yg*1e3, ez[ix_obs, :, iz_obs], '.-b')
    axc3.plot(grid.xg*1e3, ez[:, iy_obs, iz_obs], '.-r')
    axc3.set_ylim(ez_min, ez_max)
    axc3.set_xlim(-xy_max*1e3, xy_max*1e3)
    axc3.set_ylabel('$\\mathbf{E_{\\tau}\ \ [V/m]}$')
    axc3.set_xlabel('$\\mathbf{x,y\ \ [mm]}$')
    
    for ax in [axc1, axc2, axc3]:
        ax.axvline(x=grid.sigma_y_beam*1e3, color='b', linestyle='--')
        ax.axvline(x=-grid.sigma_y_beam*1e3, color='b', linestyle='--')
        ax.axvline(x=grid.sigma_x_beam*1e3, color='r', linestyle='--')
        ax.axvline(x=-grid.sigma_x_beam*1e3, color='r', linestyle='--')

    fig3.clear()
    axd1_3 = fig3.add_subplot(1,1,1)
    axd1 = plt.subplot2grid(shape=(3, 3), loc=(0,2), colspan=1, fig=fig2)

    for ax in axd1, axd1_3:
        mbl = ax.pcolormesh(grid.xg*1e3, grid.yg*1e3, np.log10(rho[:,:,iz_obs].T),
                              vmin=np.log10(rho_max)-3, vmax=np.log10(rho_max))
        ax.set_aspect('equal')
        ax.set_xlabel('$\\mathbf{x\ \ [mm]}$')
        ax.set_ylabel('$\\mathbf{y\ \ [mm]}$')
        cb = plt.colorbar(mappable=mbl, ax=ax, aspect=20)
        cb.set_label('$\\mathbf{\\log\ \\rho\ \ [e^-/m^3]}$')

    fig2.subplots_adjust(left=0.1, right=0.9,hspace=.4, wspace=.4)
    fig3.subplots_adjust(left=0.12, right=0.95,hspace=.4, wspace=.4)
    axd1_3.set_xlim(grid.xg[0]*1e3, grid.xg[-1]*1e3)
    axd1_3.set_ylim(grid.yg[0]*1e3, grid.yg[-1]*1e3)
    axd1_3.set_aspect('equal')

    axd3 = plt.subplot2grid(shape=(3, 3), loc=(2,2), colspan=1, fig=fig2)
    mbl2 = axd3.pcolormesh(grid.xg*1e3, grid.yg*1e3, phi[:,:,iz_obs].T,
                          vmin=phi_min, vmax=phi_max)
    axd3.set_aspect('equal')
    axd3.set_xlabel('$\\mathbf{x\ \ [mm]}$')
    axd3.set_ylabel('$\\mathbf{y\ \ [mm]}$')
    cb2 = plt.colorbar(mappable=mbl2, ax=axd3, aspect=20)
    cb2.set_label('$\\mathbf{\\phi\ \ [V]}$')
    
    fig2.suptitle(f'$\\tau$={grid.zg[iz_obs]:.2f} ({grid.zg[iz_obs]/grid.sigma_z_beam:.2f} $\\sigma_\\tau$)')
    fig3.suptitle(f'$\\tau$={grid.zg[iz_obs]:.2f} ({grid.zg[iz_obs]/grid.sigma_z_beam:.2f} $\\sigma_\\tau$)')
#    folder_movie='.'
    fig3.savefig(folder_movie + f'/frame_rho_{i_frame:03d}.png', dpi=150)
    for nn in 2,4,6,8:
        axd1_3.plot(nn*grid.sigma_x_beam*1e3*np.cos(t), nn*grid.sigma_y_beam*1e3*np.sin(t), color='lightgray', linestyle='-',alpha=0.3)
        axd1_3.text(nn*grid.sigma_x_beam*1e3*np.cos(np.pi/4), nn*grid.sigma_y_beam*1e3*np.sin(np.pi/4),str(nn)+'$\sigma$' , color='lightgray', alpha=0.3)
    fig2.savefig(folder_movie + f'/frame_full_{i_frame:03d}.png', dpi=150)
    fig3.savefig(folder_movie + f'/frame_rho_beam_{i_frame:03d}.png', dpi=150)

os.system(' '.join([
    'ffmpeg',
    '-i %s'%folder_movie+'/frame_rho_%03d.png',
    '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=4.*PTS"',
    '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22',
    '-codec:a aac movie_%s_rho.mp4'%tag]))

os.system(' '.join([
    'ffmpeg',
    '-i %s'%folder_movie+'/frame_full_%03d.png',
    '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=4.*PTS"',
    '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22',
    '-codec:a aac movie_%s_full.mp4'%tag]))

os.system(' '.join([
    'ffmpeg',
    '-i %s'%folder_movie+'/frame_rho_beam_%03d.png',
    '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=4.*PTS"',
    '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22',
    '-codec:a aac movie_%s_rho_beam.mp4'%tag]))

shutil.rmtree(folder_movie)
#plt.show()
