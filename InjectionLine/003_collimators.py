import numpy as np
import pysixtrack
import pickle

from cpymad.madx import Madx


mad = Madx()
mad.options.echo = False
mad.options.warn = False

mad.call('lhc_injection_fortracking.seq')
mad.use('lhcb1')

twiss = mad.twiss()
tcp1 = 'tcp.d6l7.b1'
tcp2 = 'tcp.c6l7.b1'
tcp3 = 'tcp.b6l7.b1'

nc = 5.7
epsn = 3.5e-6
gamma0 = mad.sequence['lhcb1'].beam.gamma
beta0 = mad.sequence['lhcb1'].beam.beta
epsg = epsn/beta0/gamma0

#tcp_d6l7 = pysixtrack.Elements.LimitRect(min_x = -sb_tcp1+xco_tcp1

sbx = {}
sby = {}
xco = {}
yco = {}
deg3 = 127.5
for i in range(len(twiss.name)):
    elname = twiss.name[i].split(':')[0]
    if elname in [tcp1, tcp2, tcp3]:
        sbx[elname] = np.sqrt(twiss.betx[i]*epsg)
        sby[elname] = np.sqrt(twiss.bety[i]*epsg)
        xco[elname] = twiss.x[i]
        yco[elname] = twiss.y[i]
        print(twiss.name[i])
        print(twiss.betx[i])
        print(twiss.bety[i])
        print(twiss.x[i])
        print(twiss.y[i])
        print()

tcp_d6l7 = pysixtrack.elements.LimitRect(min_x = -nc*sbx[tcp1] + xco[tcp1], max_x = nc*sbx[tcp1] + xco[tcp1])
tcp_c6l7 = pysixtrack.elements.LimitRect(min_y = -nc*sby[tcp2] + yco[tcp2], max_y = nc*sby[tcp2] + yco[tcp2])
sb3 = (sbx[tcp3]**2 + sby[tcp3]**2)**0.5
co3 = (xco[tcp3]**2 + yco[tcp3]**2)**0.5
shift1 = pysixtrack.elements.XYShift(dx=xco[tcp3],dy=yco[tcp3])
shift2 = pysixtrack.elements.XYShift(dx=-xco[tcp3],dy=-yco[tcp3])
rot1 = pysixtrack.elements.SRotation(angle=-deg3)
rot2 = pysixtrack.elements.SRotation(angle=deg3)
tcp_b6l7 = pysixtrack.elements.LimitRect(min_y = -nc*sb3, max_y = nc*sb3)



with open('line_with_ecloud_markers.pkl','rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid))

collimators = ['tcp.d6l7.b1','tcp.c6l7.b1','tcp.b6l7.b1']
##11191,5,9
for ii,elname in enumerate(line.element_names):
    if elname == tcp1:
        line.elements[ii] = tcp_d6l7

for ii,elname in enumerate(line.element_names):
    if elname == tcp2:
        line.elements[ii] = tcp_c6l7

for ii,elname in enumerate(line.element_names):
    if elname == tcp3:
        itcp3 = ii
line.elements[itcp3] = tcp_b6l7
line.insert_element(itcp3+1,shift2,'tcp.b6l7.b1.shift2')
line.insert_element(itcp3+1,rot2,'tcp.b6l7.b1.rot2')
line.insert_element(itcp3,rot1,'tcp.b6l7.b1.rot1')
line.insert_element(itcp3,shift1,'tcp.b6l7.b1.shift1')
for pp in line.element_names[itcp3-6:itcp3+6]:
    print(pp)

with open('line_with_ecloud_markers_and_collimators.pkl','wb') as fid:
    line = pickle.dump(line.to_dict(keepextra=True),fid)
