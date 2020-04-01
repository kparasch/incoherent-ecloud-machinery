from TricubicInterpolation import cTricubic as ti

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('..')

import kostas_filemanager as kfm
plt.style.use('kostas')

pinch = 'eclouds/refined_Pinch7_cut_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5'
epinch = 'eclouds/e_refined_Pinch7_cut_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5'

grid = kfm.h5_to_dict(pinch, group='grid')

Nx, Ny, Nz = grid['Nx'], grid['Ny'], grid['Nz']
dx, dy, dz = grid['dx'], grid['dy'], grid['dz']
x0, y0, z0 = grid['x0'], grid['y0'], grid['z0']
xg, yg, zg = grid['xg'], grid['yg'], grid['zg']

phi = np.empty([Nx,Ny,Nz,8])
ex = np.empty([Nx,Ny,Nz,8])
ey = np.empty([Nx,Ny,Nz,8])
ez = np.empty([Nx,Ny,Nz,8])

for kk in range(Nz):
    phi[:,:,kk,:] = kfm.h5_to_dict(pinch, group='slices/slice%d'%kk)['phi']
    ex[:,:,kk,:] = kfm.h5_to_dict(epinch, group='slices/ex_slice%d'%kk)['ex']
    ey[:,:,kk,:] = kfm.h5_to_dict(epinch, group='slices/ey_slice%d'%kk)['ey']
    ez[:,:,kk,:] = kfm.h5_to_dict(epinch, group='slices/ez_slice%d'%kk)['ez']

TI = ti.Tricubic_Interpolation(A=phi, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, method='Exact')
TIex = ti.Tricubic_Interpolation(A=ex, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, method='Exact')
TIey = ti.Tricubic_Interpolation(A=ey, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, method='Exact')
TIez = ti.Tricubic_Interpolation(A=ez, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, method='Exact')

print(max(xg),max(yg),max(zg))
xbox = max(abs(xg))
ybox = max(abs(yg))
zbox = max(abs(zg))

xmask = abs(xg) < xbox
ymask = abs(yg) < ybox
zmask = abs(zg) < zbox
xo = xg[xmask]
yo = yg[ymask]
zo = zg[zmask]

xo[-1] -= 1.e-10
yo[-1] -= 1.e-10
zo[-1] -= 1.e-10

x_obs = 2.0e-05
y_obs = 2.0e-05
z_obs = 0.00074948

xt = np.linspace(-xbox, xbox, 10000)
yt = np.linspace(-ybox, ybox, 10000)
zt = np.linspace(-zbox, zbox, 10000)

#####Ex#####
Ex_xt    = np.array([-TI.ddx(xx, y_obs, z_obs)   for xx in xt])
Ex_xo    = np.array([-TI.ddx(xx, y_obs, z_obs)   for xx in xo])
Ex_xt_ti = np.array([-TIex.val(xx, y_obs, z_obs) for xx in xt])

Ex_yt    = np.array([-TI.ddx(x_obs, yy, z_obs)   for yy in yt])
Ex_yo    = np.array([-TI.ddx(x_obs, yy, z_obs)   for yy in yo])
Ex_yt_ti = np.array([-TIex.val(x_obs, yy, z_obs) for yy in yt])

Ex_zt    = np.array([-TI.ddx(x_obs, y_obs, zz)   for zz in zt])
Ex_zo    = np.array([-TI.ddx(x_obs, y_obs, zz)   for zz in zo])
Ex_zt_ti = np.array([-TIex.val(x_obs, y_obs, zz) for zz in zt])
#####Ey#####
Ey_xt    = np.array([-TI.ddy(xx, y_obs, z_obs)   for xx in xt])
Ey_xo    = np.array([-TI.ddy(xx, y_obs, z_obs)   for xx in xo])
Ey_xt_ti = np.array([-TIey.val(xx, y_obs, z_obs) for xx in xt])

Ey_yt    = np.array([-TI.ddy(x_obs, yy, z_obs)   for yy in yt])
Ey_yo    = np.array([-TI.ddy(x_obs, yy, z_obs)   for yy in yo])
Ey_yt_ti = np.array([-TIey.val(x_obs, yy, z_obs) for yy in yt])

Ey_zt    = np.array([-TI.ddy(x_obs, y_obs, zz)   for zz in zt])
Ey_zo    = np.array([-TI.ddy(x_obs, y_obs, zz)   for zz in zo])
Ey_zt_ti = np.array([-TIey.val(x_obs, y_obs, zz) for zz in zt])
#####Ez#####
Ez_xt    = np.array([-TI.ddz(xx, y_obs, z_obs)   for xx in xt])
Ez_xo    = np.array([-TI.ddz(xx, y_obs, z_obs)   for xx in xo])
Ez_xt_ti = np.array([-TIez.val(xx, y_obs, z_obs) for xx in xt])

Ez_yt    = np.array([-TI.ddz(x_obs, yy, z_obs)   for yy in yt])
Ez_yo    = np.array([-TI.ddz(x_obs, yy, z_obs)   for yy in yo])
Ez_yt_ti = np.array([-TIez.val(x_obs, yy, z_obs) for yy in yt])

Ez_zt    = np.array([-TI.ddz(x_obs, y_obs, zz)   for zz in zt])
Ez_zo    = np.array([-TI.ddz(x_obs, y_obs, zz)   for zz in zo])
Ez_zt_ti = np.array([-TIez.val(x_obs, y_obs, zz) for zz in zt])

#####Ez#####
phi_xt    = np.array([TI.val(xx, y_obs, z_obs)   for xx in xt])
phi_xo    = np.array([TI.val(xx, y_obs, z_obs)   for xx in xo])

phi_yt    = np.array([TI.val(x_obs, yy, z_obs)   for yy in yt])
phi_yo    = np.array([TI.val(x_obs, yy, z_obs)   for yy in yo])

phi_zt    = np.array([TI.val(x_obs, y_obs, zz)   for zz in zt])
phi_zo    = np.array([TI.val(x_obs, y_obs, zz)   for zz in zo])




maxEx = 1.1*max([max(abs(Ex_xt)),max(abs(Ex_yt)),max(abs(Ex_zt))])
fig = plt.figure(1,figsize=[18,5])
ax1 = fig.add_subplot(131)
ax1.plot(xt, Ex_xt,'r')
ax1.plot(xo, Ex_xo,'b.')
ax1.plot(xt, Ex_xt_ti,'k--')
ax1.set_xlabel('x')
ax1.set_ylabel('Ex')
ax1.set_ylim(-maxEx, maxEx)

ax2 = fig.add_subplot(132)
ax2.plot(yt, Ex_yt,'r')
ax2.plot(yo, Ex_yo,'b.')
ax2.plot(yt, Ex_yt_ti,'k--')
ax2.set_xlabel('y')
ax2.set_ylim(-maxEx, maxEx)

ax3 = fig.add_subplot(133)
ax3.plot(zt, Ex_zt,'r')
ax3.plot(zo, Ex_zo,'b.')
ax3.plot(zt, Ex_zt_ti,'k--')
ax3.set_xlabel('z')
ax3.set_ylim(-maxEx, maxEx)

maxEy = 1.1*max([max(abs(Ey_xt)),max(abs(Ey_yt)),max(abs(Ey_zt))])
fig2 = plt.figure(2,figsize=[18,5])
ax21 = fig2.add_subplot(131)
ax21.plot(xt, Ey_xt,'r')
ax21.plot(xo, Ey_xo,'b.')
ax21.plot(xt, Ey_xt_ti,'k--')
ax21.set_xlabel('x')
ax21.set_ylabel('Ey')
ax21.set_ylim(-maxEy, maxEy)

ax22 = fig2.add_subplot(132)
ax22.plot(yt, Ey_yt,'r')
ax22.plot(yo, Ey_yo,'b.')
ax22.plot(yt, Ey_yt_ti,'k--')
ax22.set_xlabel('y')
ax22.set_ylim(-maxEy, maxEy)

ax23 = fig2.add_subplot(133)
ax23.plot(zt, Ey_zt,'r')
ax23.plot(zo, Ey_zo,'b.')
ax23.plot(zt, Ey_zt_ti,'k--')
ax23.set_xlabel('z')
ax23.set_ylim(-maxEy, maxEy)


maxEz = 1.1*max([max(abs(Ez_xt)),max(abs(Ez_yt)),max(abs(Ez_zt))])
fig3 = plt.figure(3,figsize=[18,5])
ax31 = fig3.add_subplot(131)
ax31.plot(xt, Ez_xt,'r')
ax31.plot(xo, Ez_xo,'b.')
ax31.plot(xt, Ez_xt_ti,'k--')
ax31.set_xlabel('x')
ax31.set_ylabel('Ez')
ax31.set_ylim(-maxEz, maxEz)

ax32 = fig3.add_subplot(132)
ax32.plot(yt, Ez_yt,'r')
ax32.plot(yo, Ez_yo,'b.')
ax32.plot(yt, Ez_yt_ti,'k--')
ax32.set_xlabel('y')
ax32.set_ylim(-maxEz, maxEz)

ax33 = fig3.add_subplot(133)
ax33.plot(zt, Ez_zt,'r')
ax33.plot(zo, Ez_zo,'b.')
ax33.plot(zt, Ez_zt_ti,'k--')
ax33.set_xlabel('z')
ax33.set_ylim(-maxEz, maxEz)

maxPhi = 1.1*max([max(abs(phi_xt)),max(abs(phi_yt)),max(abs(phi_zt))])
fig4 = plt.figure(4,figsize=[18,5])
ax41 = fig4.add_subplot(131)
ax41.plot(xt, phi_xt,'r')
ax41.plot(xo, phi_xo,'b.')
ax41.set_xlabel('x')
ax41.set_ylabel('phi')
ax41.set_ylim(-maxPhi, maxPhi)

ax42 = fig4.add_subplot(132)
ax42.plot(yt, phi_yt,'r')
ax42.plot(yo, phi_yo,'b.')
ax42.set_xlabel('y')
ax42.set_ylim(-maxPhi, maxPhi)

ax43 = fig4.add_subplot(133)
ax43.plot(zt, phi_zt,'r')
ax43.plot(zo, phi_zo,'b.')
ax43.set_xlabel('z')
ax43.set_ylim(-maxPhi, maxPhi)
plt.show()




