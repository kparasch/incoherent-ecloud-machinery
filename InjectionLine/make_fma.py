import sys 
import matplotlib.cm as cm
sys.path.append('../Tools')
import kostas_filemanager as kfm 
import matplotlib.pyplot as plt 
import numpy as np
import tune_diagram

#AA = kfm.h5_to_obj('fma_tunes_ec_scale0.10.h5')
#AA = kfm.h5_to_obj('fma_tunes_ec_scale0.05_IMO_0.h5')
AA = kfm.h5_to_obj('fma_tunes_ec_scale0.10_IMO_0_ref.h5')
#AA = kfm.h5_to_obj('fma_tunes_ec_scale0.05.h5')


#plt.plot(AA.qx[0,:], AA.qy[0,:],'.')
#plt.plot(AA.q1[0,:], AA.q2[0,:],'.')

last = -2
mq1 = np.mean(AA.q1[:last,:], axis=0)
mq2 = np.mean(AA.q2[:last,:], axis=0)
q1_std = np.std(AA.q1[:last,:], axis=0)
q2_std = np.std(AA.q2[:last,:], axis=0)

dq = np.sqrt(q1_std**2 + q2_std**2)
dq = np.sqrt((AA.q1[-1,:] - AA.q1[-2,:])**2 + (AA.q2[-1,:] - AA.q2[-2,:])**2)

vmin = -7
vmax = -3

qx_range = [0.2,0.4]
qy_range = [0.2,0.4]
qx_range = [0.268,0.284]
qy_range = [0.294,0.3125]
orders = 40
plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
#RL = tune_diagram.ResonanceLines([0.25,0.3], [0.27,0.33], range(orders), 1)
RL = tune_diagram.ResonanceLines(qx_range, qy_range, range(1,orders+1), 1)
RL.plot_resonance(fig1)
cb = ax1.scatter(mq1, mq2,c=np.log10(dq), s=5, vmin=vmin, vmax=vmax, cmap=cm.jet)
fig1.colorbar(cb)
plt.figure(2)
plt.scatter(np.sqrt(2*AA.J1), np.sqrt(2*AA.J2),c=np.log10(dq), s=5, vmin=vmin, vmax=vmax, cmap=cm.jet)
plt.colorbar()
plt.show()

