import pysixtrack
import matplotlib.pyplot as plt
import os
import shutil
import sys
sys.path.append('../Tools/')
from replaceline import replaceline_and_save

from cpymad.madx import Madx

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--noblock', dest='plot_block', action='store_false')
parser.add_argument('--qx0', nargs='?', default=62.270, type=float)
parser.add_argument('--qy0', nargs='?', default=60.295, type=float)
parser.add_argument('--qprime', nargs='?', default=20., type=float)
parser.add_argument('--I_MO', nargs='?', default=40., type=float)
parser.add_argument('--VRF400', nargs='?', default=8., type=float)
args = parser.parse_args()

mad = Madx()


fname = 'lhc2018_injection.mask'
replaceline_and_save(fname,
                     findln='qx0 =',
                     newline=f'qx0 = {args.qx0:.3f};'
                    )
replaceline_and_save(fname,
                     findln='qy0 =',
                     newline=f'qy0 = {args.qy0:.3f};'
                    )
replaceline_and_save(fname,
                     findln='qprime =',
                     newline=f'qprime = {args.qprime:.1f};'
                    )
replaceline_and_save(fname,
                     findln='I_MO =',
                     newline=f'I_MO = {args.I_MO:.1f};'
                    )
replaceline_and_save(fname,
                     findln='VRF400:=',
                     newline=f'VRF400:={args.VRF400:.1f};'
                    )

mad.call('lhc2018_injection.mask')

#mad.command.esave(file='lattice_errors.err')

os.remove('db4')
os.remove('db5')
os.remove('fidel')
os.remove('slhc')
os.remove('wise')
os.remove('twiss.b1')
os.remove('twiss.b2')
shutil.rmtree('temp')

mad.use('lhcb1')
twiss = mad.twiss(sequence='lhcb1')

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(twiss.s, twiss.x)
ax1.set_xlabel('$\mathbf{s [m]}$')
ax1.set_ylabel('$\mathbf{x_{CO} [m]}$')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)

ax2.plot(twiss.s, twiss.y)
ax2.set_xlabel('$\mathbf{s [m]}$')
ax2.set_ylabel('$\mathbf{y_{CO} [m]}$')

plt.show(block=args.plot_block)

