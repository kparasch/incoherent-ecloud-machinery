import pysixtrack
import matplotlib.pyplot as plt
import os
import shutil

from cpymad.madx import Madx

mad = Madx()

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

plt.show()

