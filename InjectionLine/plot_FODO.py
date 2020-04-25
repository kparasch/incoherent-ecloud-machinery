from cpymad.madx import Madx
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('kostas')

plt.close('all')

mad = Madx()
mad.options.echo = False
mad.options.warn = False

mad.call('lhc_injection_fortracking.seq')
mad.use('lhcb1')

twiss = mad.twiss()

draw_ec = True

s_cell = 1240
e_cell = 1400
mq0 = 1251 #defocusing
mq1 = 1293 #focusing
mq2 = 1341 #defocusing
mq3 = 1383 #focusing
mq_list = [1251, 1293, 1341, 1383]
mb_list = [1219, 1227, 1239, 1261, 1273, 1281, 1309, 1317, 1329, 1351, 1363, 1371, 1399]
mq_HL = 3.1/2.
mb_HL = 14.3/2.

mq_h=0.1
mb_h=0.05
fig1 = plt.figure(1, figsize=[9,5])
ax1 = fig1.add_subplot(111)
[ax1.axvspan(twiss.s[mq] - mq_HL, twiss.s[mq] + mq_HL, ymax=mq_h, facecolor='r', alpha=0.5) for mq in mq_list]
[ax1.axvspan(twiss.s[mb] - mb_HL, twiss.s[mb] + mb_HL, ymax=mb_h, facecolor='b', alpha=0.5) for mb in mb_list]
ax1.axvspan(twiss.s[mb_list[0]] - mb_HL, twiss.s[mb_list[0]] + mb_HL, ymax=0, facecolor='b', alpha=0.5, label='$\\mathbf{MB}$') 
ax1.axvspan(twiss.s[mq_list[0]] - mq_HL, twiss.s[mq_list[0]] + mq_HL, ymax=0, facecolor='r', alpha=0.5, label='$\\mathbf{MQ}$') 
ax1.plot(twiss.s[s_cell:e_cell+1], twiss.betx[s_cell:e_cell+1],'k-', label='$\\mathbf{\\beta_x\ [m]}$')
ax1.plot(twiss.s[s_cell:e_cell+1], twiss.bety[s_cell:e_cell+1],'r-', label='$\\mathbf{\\beta_y\ [m]}$')
ax1.plot([twiss.s[s_cell], twiss.s[e_cell]],[-10,-10],'b-', label='$\\mathbf{D_x\ [m]}$')
ax2 = ax1.twinx()
ax2.plot(twiss.s[s_cell:e_cell+1], twiss.dx[s_cell:e_cell+1],'b-')
ax2.set_ylim(0,3.0)
ax1.set_ylim(0, 200)
ax1.set_xlabel('$\\mathbf{s\ [m]}$')
ax1.set_ylabel('$\\mathbf{\\beta\ [m]}$')
ax2.set_ylabel('$\\mathbf{D\ [m]}$')
if draw_ec:
    ax1.arrow(twiss.s[1363], y=130, dx=0, dy=-20, width=1.1, color='g')
    ax1.arrow(twiss.s[1317], y=130, dx=0, dy=-20, width=1.1, color='g')
    ax1.arrow(twiss.s[1273], y=130, dx=0, dy=-20, width=1.1, color='g')
    ax1.axvline(twiss.s[1363], color='g', linestyle='dashed')
    ax1.axvline(twiss.s[1317], color='g', linestyle='dashed')
    ax1.axvline(twiss.s[1273], color='g', linestyle='dashed')
ax1.legend(loc='upper right', fontsize=20, bbox_to_anchor=(1.0, 1.00), bbox_transform=plt.gcf().transFigure)
ax1.set_xlim(twiss.s[s_cell], twiss.s[e_cell])
ax1.grid(False)
ax2.grid(False)
fig1.subplots_adjust(right=0.69)
ax1.ticklabel_format(style='plain')
ax2.ticklabel_format(style='plain')
#[plt.axvline(twiss.s[mq] - mq_HL) for mq in mq_list]
#[plt.axvline(twiss.s[mq] + mq_HL) for mq in mq_list]
#[plt.axvline(twiss.s[mb] - mb_HL) for mb in mb_list]
#[plt.axvline(twiss.s[mb] + mb_HL) for mb in mb_list]
plt.show()

