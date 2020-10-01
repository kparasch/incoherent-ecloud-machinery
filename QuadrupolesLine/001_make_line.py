import pickle
import pysixtrack
import numpy as np
from cpymad.madx import Madx
import matplotlib.pyplot as plt

import sys
sys.path.append('../Tools')
import pyht_beamsize
#plt.style.use('kostas')
plt.close('all')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--noblock', dest='plot_block', action='store_false')
args = parser.parse_args()

mad = Madx()

mad.options.echo = False
mad.options.warn = False
#mad.options.info = False

mad.call('lhc_injection_fortracking.seq')
mad.use('lhcb1')


twiss = mad.twiss()

option = 'every-mb.b'
#option = 'between-mq'

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(twiss.s, twiss.x)
ax1.set_xlabel('$\mathbf{s [m]}$')
ax1.set_ylabel('$\mathbf{x_{CO},y_{CO} [m]}$')
ax1.plot(twiss.s, twiss.y)
ax1.set_xlabel('$\mathbf{s [m]}$')
ax1.set_ylabel('$\mathbf{y_{CO} [m]}$')

#find dipoles and quadrupoles between  s.arc and e.arc
name_arcs = ['34','45','56','67','78','81','12','23']
#Find MBB in arcs
start_arcs = []
end_arcs = []
for i,el in enumerate(twiss.name):
    if el[0:5] == 's.arc':
        start_arcs.append(i)
    if el[0:5] == 'e.arc':
        end_arcs.append(i)

mbb_arc = []
for i in range(8): #arcs
    mbb_list = []
    for j in range(start_arcs[i],end_arcs[i]):
        el_name = twiss.name[j].split(':')[0]
        if len(el_name.split('..')) > 1:
            continue
        if el_name[0:4] in ['mb.b']:
            mbb_list.append(j)
    mbb_arc.append(mbb_list)

#Find MQF and MQD
#start_names = ['ip3:1', 'ip4:1', 'ip5:1', 'ip6:1', 'ip7:1', 'ip8:1', 'ip1:1', 'ip2:1']
#end_names   = ['ip4:1', 'ip5:1', 'ip6:1', 'ip7:1', 'ip8:1', 'ip1:1', 'ip2:1','lhcb1$end:1']
#start_arcs = []
#end_arcs = []
#jjs=0
#jje=0
#for i,el in enumerate(twiss.name):
#    if jjs < 8 and el == start_names[jjs]:
#        start_arcs.append(i)
#        jjs+=1
#    if jje < 8 and el == end_names[jje]:
#        end_arcs.append(i)
#        jje+=1

mqf_arc = []
mqd_arc = []
for i in range(8): #arcs
    mqf_list = []
    mqd_list = []
    for j in range(start_arcs[i],end_arcs[i]):
        el_name = twiss.name[j].split(':')[0]
        if len(el_name.split('..')) > 1:
            continue
        if el_name.split('.')[0] in ['mq']:
            assert el_name + '..1' == twiss.name[j-2].split(':')[0]
            focusing_quad = twiss.k1l[j-2] > 0
            if focusing_quad:
                mqf_list.append(j)
            else:
                mqd_list.append(j)
    mqf_arc.append(mqf_list)
    mqd_arc.append(mqd_list)

print([len(i) for i in mqf_arc])
print([len(i) for i in mqd_arc])

#for ii in range(8):
#    print(f'Arc: {ii}')
#    for mqf in mqf_arc[ii]:
#        print(twiss.name[mqf])
#    for mqd in mqd_arc[ii]:
#        print(twiss.name[mqd])


ecl = {'mbb' : 14.3*3.,
       'mqf' : 3.1,
       'mqd' : 3.1}

arc = {'mbb' : mbb_arc,
       'mqf' : mqf_arc,
       'mqd' : mqd_arc}

ecloud_arc = {}
ecloud_lengths = {}
for ec_type in ['mbb', 'mqf', 'mqd']:
    ecloud_arc[ec_type] = []
    ecloud_lengths[ec_type] = []
    for ec_list in arc[ec_type]:
        ec_s = np.array(twiss.s[ec_list])
        ecloud_s = ec_s.copy()
        ecloud_l = ecl[ec_type]*np.ones_like(ecloud_s)
        ecloud_arc[ec_type].append(list(ecloud_s))
        ecloud_lengths[ec_type].append(list(ecloud_l))

mad.input('seqedit, sequence=lhcb1;flatten;')
ecloud_lengths_dict = {}
ec_types = ['mbb', 'mqf','mqd']
for ec_type in ec_types:
    print('Number of eclouds per arc: %d'%len(ecloud_arc[ec_type][0]))
    for i,eclouds in enumerate(ecloud_arc[ec_type]):
        ecloud_l = ecloud_lengths[ec_type][i]
        for j,ec in enumerate(eclouds):
            ecloud_length = ecloud_l[j]
            ecloud_name = f'ecloud.{ec_type}.{name_arcs[i]}.{j}'
            ecloud_lengths_dict[ecloud_name] = ecloud_length
            mad.input((f'install, element={ecloud_name}, class=marker, at={ec};'))
            #mad.input(f'install, ecloud.{i%d}.{j%d}
mad.input('flatten;endedit;')

mad.use('lhcb1')

line = pysixtrack.Line.from_madx_sequence(mad.sequence['lhcb1'])
#for el,name in zip(line.elements, line.element_names):
#    if name[0:6] == 'ecloud':
#        el._extra = ['ecloud_length'] # : ecloud_lengths_dict[name]}]
#        #el._fields.append('ecloud_length')# : ecloud_lengths_dict[name]}]
#        setattr(el,'ecloud_length', ecloud_lengths_dict[name])
#        #el._extra_fields = [{'ecloud_length' : ecloud_lengths_dict[name]}]
#line._extra.append(ecloud_lengths_dict)

with open('line_with_ecloud_markers.pkl','wb') as fid:
    pickle.dump(line.to_dict(keepextra=True), fid)

with open('ecloud_lengths.pkl','wb') as fid:
    pickle.dump(ecloud_lengths_dict, fid)

seq = 'lhcb1'
mad.use(seq)
twiss_table = mad.twiss()
optics = {'betx'      : twiss_table.betx[0],
          'bety'      : twiss_table.bety[0],
          'alfx'      : twiss_table.alfx[0],
          'alfy'      : twiss_table.alfy[0],
          'q1'        : mad.table['summ']['q1'][0],
          'q2'        : mad.table['summ']['q2'][0],
          'dq1'       : mad.table['summ']['dq1'][0],
          'dq2'       : mad.table['summ']['dq2'][0],
          'rf_volt_V'   : sum(twiss_table.volt)*1.e6,
          'rf_freq_Hz'   : max(twiss_table.freq)*1.e6,
          'rf_harmon' : max(twiss_table.harmon),
          'rf_lag'    : twiss_table.lag[twiss_table.harmon != 0][0],
          'length'    : twiss_table.s[-1],
          'alfa'      : twiss_table.alfa[0], 
          'beta0'     : mad.sequence[seq].beam.beta,
          'gamma0'    : mad.sequence[seq].beam.gamma,
          'p0c_eV'    : mad.sequence[seq].beam.pc*1.e9,
         }

for kk in optics:
    print(f'{kk} = {optics[kk]}')

with open('optics.pkl','wb') as fid:
    pickle.dump(optics, fid)
          #'M'         : M,
          #'Ms'        : Ms,
          #'W'         : W,
          #'invW'      : invW,
          #'R'         : R

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
mad.input('seqedit, sequence=lhcb1;flatten;')
mad.input('cycle,start=ip5;')
mad.input('flatten;endedit;')
mad.use('lhcb1')
survey = mad.survey()
twiss = mad.twiss()
ax3.plot(survey.z, survey.x,'b')
ax3.set_axis_off()
ax3.set_aspect('equal')

center_z = np.mean(survey.z)
center_x = np.mean(survey.x)
roff = 0.5e3
betx_ec = []
bety_ec = []
dispx_ec = []
dispy_ec = []
name_ec = []
name_ip = []
s_ec = []
sig11_ec = []
sig33_ec = []
s_ip = []

epsn_x = mad.sequence['lhcb1'].beam.exn
epsn_y = mad.sequence['lhcb1'].beam.eyn
sigma_z = mad.sequence['lhcb1'].beam.sigt

for i in range(len(survey.name)):
    if survey.name[i][0:3] == 'ecl':
        ax3.plot(survey.z[i], survey.x[i], 'r.')
        betx_ec.append(twiss.betx[i])
        bety_ec.append(twiss.bety[i])
        dispx_ec.append(twiss.dx[i])
        dispy_ec.append(twiss.dy[i])
        name_ec.append(twiss.name[i])
        s_ec.append(twiss.s[i])
        sig11_ec.append(twiss.sig11[i])
        sig33_ec.append(twiss.sig33[i])

    if survey.name[i][0:2] == 'ip' and survey.name[i][3] == ':':

        s_ip.append(survey.s[i])
        name_ip.append(survey.name[i][0:3].upper())
        ax3.plot(survey.z[i], survey.x[i], 'ko')

        text_z = survey.z[i] - center_z
        text_x = survey.x[i] - center_x
        text_r = np.sqrt((text_x)**2+(text_z)**2) - roff
        text_theta = np.arctan2(text_x,text_z)
        text_z = text_r* np.cos(text_theta) + center_z
        text_x = text_r* np.sin(text_theta) + center_x
        ax3.text(text_z, text_x, '$\\mathbf{' + survey.name[i][0:3].upper() + '}$', 
                 horizontalalignment='center', verticalalignment='center')
        print(survey.name[i])

latex_name_ip = ['$\\mathbf{'+name+'}$' for name in name_ip]
sig1_nonlin, sig2_nonlin, sig1_lin, sig2_lin = pyht_beamsize.sigmas(epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z, 
          beta_x=np.mean(betx_ec), beta_y=np.mean(bety_ec), D_x=np.mean(dispx_ec), D_y=np.mean(dispy_ec),
          alpha_mom_compaction=optics['alfa'], circumference=optics['length'],
          rf_harmonic=optics['rf_harmon'], V_rf=optics['rf_volt_V'], gamma0=optics['gamma0'])


sx = np.mean(np.sqrt(sig11_ec))*1000
sxerr = np.std(np.sqrt(sig11_ec))*1000
sy = np.mean(np.sqrt(sig33_ec))*1000
syerr = np.std(np.sqrt(sig33_ec))*1000
print(f'Horizontal beam size on eclouds = {sx:.4f} +- {sxerr:.4f} ({100*sxerr/sx:.1f}%) [mm]')
print(f'Vertical beam size on eclouds = {sy:.4f} +- {syerr:.4f} ({100*syerr/sy:.1f}%) [mm]')
print(f'Max horizontal beam size on eclouds = {1000*np.max(np.sqrt(sig11_ec))} [mm]')
print(f'Max vertical beam size on eclouds = {1000*np.max(np.sqrt(sig33_ec))} [mm]')
print('PyHEADTAIL beam size:')
print(f'  Linear bucket:')
print(f'    Horizontal beam size = {sig1_lin*1000:.4f} [mm]')
print(f'    Vertical   beam size = {sig2_lin*1000:.4f} [mm]')
print(f'  Non-linear bucket:')
print(f'    Horizontal beam size = {sig1_nonlin*1000:.4f} [mm]')
print(f'    Vertical   beam size = {sig2_nonlin*1000:.4f} [mm]')
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)

ax4.plot(s_ec,betx_ec,'k.-',label='$\\mathbf{\\beta_x}$')
ax4.plot(s_ec,bety_ec,'r.-',label='$\\mathbf{\\beta_y}$')
ax4.plot([0,1],[0,0],'b.-', label='$\\mathbf{D_x}$')
ax42 = ax4.twinx()
ax42.plot(s_ec,dispx_ec,'b.-',label='$\\mathbf{D_x}$')
[ax4.axvline(this_s_ip,c='k', linestyle='dashed', alpha=0.5) for this_s_ip in s_ip]
ax4.set_ylabel('$\\mathbf{\\beta\ [m]}$')
ax42.set_ylabel('$\\mathbf{D\ [m]}$',color='b')
ax4.set_ylim(0,200.)
ax42.set_ylim(1,3.5)
ax4.set_xlabel('$\\mathbf{s\ [m]}$')
plt.xticks(s_ip, latex_name_ip)
ax4.legend(loc='upper left')
ax4.grid(False)
ax42.grid(False)
fig4.subplots_adjust(right=0.90)

fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111)
ax5.plot(s_ec,np.sqrt(sig11_ec),'k.-',label='$\\mathbf{\\sigma_x}$')
ax5.plot(s_ec,np.sqrt(sig33_ec),'r.-',label='$\\mathbf{\\sigma_y}$')
ax5.axhline(sig1_lin, label='$\\mathbf{PyHT\ linear}$', color='b', linewidth=3.0)
ax5.axhline(sig2_lin, color='b', linewidth=3.0)
ax5.axhline(sig1_nonlin, label='$\\mathbf{PyHT\ nonlinear}$', color='g', linewidth=3.0)
#ax5.axhline(sig1_nonlin, label='$\\mathbf{PyHT\ non\\mbox{-}linear}$', color='g', linewidth=3.0)
ax5.axhline(sig2_nonlin, color='g', linewidth=3.0)
[ax5.axvline(this_s_ip,c='k', linestyle='dashed', alpha=0.5) for this_s_ip in s_ip]
ax5.legend()
ax5.set_ylabel('$\\mathbf{\\sigma\ [m]}$')
ax5.set_xlabel('$\\mathbf{s\ [m]}$')
ax5.grid(False)
plt.xticks(s_ip, latex_name_ip)
plt.show(block=args.plot_block)
