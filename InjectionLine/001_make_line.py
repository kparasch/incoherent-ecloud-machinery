import pickle
import pysixtrack
import numpy as np
from cpymad.madx import Madx
import matplotlib.pyplot as plt

plt.style.use('kostas')
plt.close('all')

mad = Madx()

mad.call('lhc_injection_fortracking.seq')
mad.use('lhcb1')


twiss = mad.twiss()

#fig1 = plt.figure(1)
#ax1 = fig1.add_subplot(111)
#ax1.plot(twiss.s, twiss.x)
#ax1.set_xlabel('$\mathbf{s [m]}$')
#ax1.set_ylabel('$\mathbf{x_{CO},y_{CO} [m]}$')
#ax1.plot(twiss.s, twiss.y)
#ax1.set_xlabel('$\mathbf{s [m]}$')
#ax1.set_ylabel('$\mathbf{y_{CO} [m]}$')
#
#fig2 = plt.figure(2)
#ax2 = fig2.add_subplot(111)
#
#ax2.plot(twiss.s, twiss.betx)
#ax2.plot(twiss.s, twiss.bety)
#ax2.plot(twiss.s, twiss.dx)
#ax2.set_xlabel('$\mathbf{s [m]}$')
#ax2.set_ylabel('$\mathbf{\\beta_{x},\\beta_{y},D_{x} [m]}$')

#find quadrupoles between  s.arc and e.arc
start_arcs = []
end_arcs = []
name_arcs = []
for i,el in enumerate(twiss.name):
    if el[0:5] == 's.arc':
        start_arcs.append(i)
        name_arcs.append(el[6:8])
    if el[0:5] == 'e.arc':
        end_arcs.append(i)

mq_arc = []
for i in range(8): #arcs
    mq_list = []
    for j in range(start_arcs[i],end_arcs[i]):
        el_name = twiss.name[j].split(':')[0]
        if len(el_name.split('..')) > 1:
            continue
        if el_name.split('.')[0] in ['mq']:
            mq_list.append(j)
    mq_arc.append(mq_list)

print([len(i) for i in mq_arc])

ecloud_arc = []
ecloud_lengths = []
for mq_list in mq_arc:
    mq_s = np.array(twiss.s[mq_list])
    ecloud_s = (mq_s[:-1] + mq_s[1:])/2.
    ecloud_l = (mq_s[1:] - mq_s[:-1])
    ecloud_arc.append(list(ecloud_s))
    ecloud_lengths.append(list(ecloud_l))

print('Number of eclouds per arc: %d'%len(ecloud_arc[0]))
mad.input('seqedit, sequence=lhcb1;flatten;')
for i,eclouds in enumerate(ecloud_arc):
    for j,ec in enumerate(eclouds):
        ecloud_name = f'ecloud.{name_arcs[i]}.{j}'
        mad.input((f'install, element={ecloud_name}, class=marker, at={ec};'))
        #mad.input(f'install, ecloud.{i%d}.{j%d}
mad.input('flatten;endedit;')

mad.use('lhcb1')

line = pysixtrack.Line.from_madx_sequence(mad.sequence['lhcb1'])
with open('line_with_ecloud_markers.pkl','wb') as fid:
    pickle.dump(line.to_dict(keepextra=True), fid)

#mq_s = []
#for i in range(len(twiss.name)):
#    el_name = twiss.name[i].split(':')[0]
#    split_name = el_name.split('.')
#
#    el_magnet = split_name[0]
#    if len(split_name) > 1:
#        el_cell = split_name[1]
#    if len(split_name) > 2:
#        el_beam = split_name[2]
#    el_s = twiss.s[i]
#    if el_name[0] in ['s','e']:
#        print(el_name, el_s, i)
#    if i > 20624 and i < 22724:
#        if len(el_name.split('..')) > 1: 
#            continue
#        if el_magnet in ['mq','mb']:
#            if el_magnet == 'mq':
#                mq_s.append(el_s)
#            print(el_name, el_s)

#mq_s = np.array(mq_s)
#print(mq_s[1:] - mq_s[:-1])
#print((mq_s[1:] - mq_s[:-1])/2. + mq_s[:-1])

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
name_ec = []
name_ip = []
s_ec = []
sig11_ec = []
sig33_ec = []
s_ip = []
for i in range(len(survey.name)):
    if survey.name[i][0:3] == 'ecl':
        ax3.plot(survey.z[i], survey.x[i], 'r.')
        betx_ec.append(twiss.betx[i])
        bety_ec.append(twiss.bety[i])
        dispx_ec.append(twiss.dx[i])
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


fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)

ax4.plot(s_ec,betx_ec,'.-')
ax4.plot(s_ec,bety_ec,'.-')
ax4.plot(s_ec,dispx_ec,'.-')
[ax4.axvline(this_s_ip,c='k') for this_s_ip in s_ip]
plt.xticks(s_ip, name_ip)

fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111)
ax5.plot(s_ec,np.sqrt(sig11_ec),'.-')
ax5.plot(s_ec,np.sqrt(sig33_ec),'.-')
[ax5.axvline(this_s_ip,c='k') for this_s_ip in s_ip]
plt.xticks(s_ip, name_ip)
plt.show()
