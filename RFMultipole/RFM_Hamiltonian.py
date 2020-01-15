import numpy as np
from scipy.constants import c as clight
import sympy
import math
from TricubicInterpolation import cTricubic as cTI
import matplotlib.pyplot as plt


class RFMultipole:
    knl = [0]
    ksl = [0]
    pn = [0]
    ps = [0]
    frequency = 0.0
    voltage = 0.0
    lag = 0.0

    def init_H(self):
        self.Hamiltonian = 0.0
        order = max(len(self.knl), len(self.ksl)) - 1
        while len(self.knl) < order + 1:
            self.knl.append(0)
        while len(self.ksl) < order + 1:
            self.ksl.append(0)
        while len(self.pn) < order + 1:
            self.pn.append(0)
        while len(self.ps) < order + 1:
            self.ps.append(0)
        k = 2.0 * np.pi * self.frequency / clight
        x = sympy.Symbol("x", real=True)
        y = sympy.Symbol("y", real=True)
        tau = sympy.Symbol("tau", real=True)
        self.xsymb = x
        self.ysymb = y
        self.tausymb = tau
        ktau = k * tau
        deg2rad = np.pi / 180.0
        self.pn = [deg2rad * pp for pp in self.pn]
        self.ps = [deg2rad * pp for pp in self.ps]

        for ii in range(order + 1):
            pn_ii = self.pn[ii] - ktau
            ps_ii = self.ps[ii] - ktau
            cn = sympy.cos(pn_ii)
            sn = sympy.sin(pn_ii)
            cs = sympy.cos(ps_ii)
            ss = sympy.sin(ps_ii)
            expr = self.knl[ii] * cn + sympy.I * self.ksl[ii] * cs
            expr = expr * ((x + sympy.I * y) ** (ii + 1))
            expr = sympy.re(expr)
            expr /= 1.0 * math.factorial(ii + 1)
            self.Hamiltonian += expr
        self.Hamiltonian = sympy.simplify(self.Hamiltonian)
        self.d_Hamiltonian_dx = self.Hamiltonian.diff(x)
        self.d_Hamiltonian_dy = self.Hamiltonian.diff(y)
        self.d_Hamiltonian_dtau = self.Hamiltonian.diff(tau)
        self.d_Hamiltonian_dxdy = self.Hamiltonian.diff(x, y)
        self.d_Hamiltonian_dxdtau = self.Hamiltonian.diff(x, tau)
        self.d_Hamiltonian_dydtau = self.Hamiltonian.diff(y, tau)
        self.d_Hamiltonian_dxdydtau = self.Hamiltonian.diff(x, y, tau)

        self.H = sympy.lambdify((x,y,tau), self.Hamiltonian, modules='numpy')
        self.dHdx = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dx, modules='numpy')
        self.dHdy = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dy, modules='numpy')
        self.dHdtau = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dtau, modules='numpy')
        self.dHdxdy = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dxdy, modules='numpy')
        self.dHdxdtau = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dxdtau, modules='numpy')
        self.dHdydtau = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dydtau, modules='numpy')
        self.dHdxdydtau = sympy.lambdify((x,y,tau), self.d_Hamiltonian_dxdydtau, modules='numpy')
        print('Hamiltonian is:')
        print(self.Hamiltonian)

    def get_A(self, xmax, ymax, zmax, nx, ny, nz):
        
        xg = np.linspace(-xmax, xmax, nx)
        yg = np.linspace(-ymax, ymax, ny)
        zg = np.linspace(-zmax, zmax, nz)

        A = np.empty([nx,ny,nz,8])
        for ii, xx in enumerate(xg): 
            for jj, yy in enumerate(yg): 
                A[ii, jj, :, 0] = self.H(xx, yy, zg)
                A[ii, jj, :, 1] = self.dHdx(xx, yy, zg)
                A[ii, jj, :, 2] = self.dHdy(xx, yy, zg)
                A[ii, jj, :, 3] = self.dHdtau(xx, yy, zg)
                A[ii, jj, :, 4] = self.dHdxdy(xx, yy, zg)
                A[ii, jj, :, 5] = self.dHdxdtau(xx, yy, zg)
                A[ii, jj, :, 6] = self.dHdydtau(xx, yy, zg)
                A[ii, jj, :, 7] = self.dHdxdydtau(xx, yy, zg)

        x0 = xg[0]
        y0 = yg[0]
        z0 = zg[0]
        dx = xg[1] - xg[0]
        dy = yg[1] - yg[0]
        dz = zg[1] - zg[0]

        return A, x0, y0, z0, dx, dy, dz

    def report(self, xmax, ymax, zmax, nx, ny, nz):

        A, x0, y0, z0, dx, dy, dz = self.get_A( xmax, ymax, zmax, nx, ny, nz )

        TI = cTI.Tricubic_Interpolation(
            A=A, dx=dx, dy=dy, dz=dz, x0=x0, y0=y0, z0=z0, method="Exact"
        )

        xobs = x0 + dx * (nx//2 + 2)
        yobs = y0 + dy * (ny//2 + 2)
        zobs = z0 + dz * (nz//2 + 2)

        xg = np.linspace(-xmax, xmax, nx)[:-1]
        yg = np.linspace(-ymax, ymax, ny)[:-1]
        zg = np.linspace(-zmax, zmax, nz)[:-1]
        xl = np.linspace(-xmax, xmax, 1000)[:-1]
        yl = np.linspace(-ymax, ymax, 1000)[:-1]
        zl = np.linspace(-zmax, zmax, 1000)[:-1]

        xobs = xg[nx//2 + 2]
        yobs = yg[ny//2 + 2]
        zobs = zg[nz//2 + 2]

        fig = plt.figure(figsize=[18,18])
        ax11 = fig.add_subplot(3,3,1)
        ax11.plot(xl, self.dHdx(xl, yobs, zobs), 'b-')
        ax11.plot(xl, np.array([ TI.ddx(xx, yobs, zobs) for xx in xl]),'r--')
        ax11.plot(xg, np.array([ TI.ddx(xx, yobs, zobs) for xx in xg]),'ro')

        ax12 = fig.add_subplot(3,3,2)
        ax12.plot(yl, self.dHdx(xobs, yl, zobs), 'b-')
        ax12.plot(yl, np.array([ TI.ddx(xobs, yy, zobs) for yy in yl]),'r--')
        ax12.plot(yg, np.array([ TI.ddx(xobs, yy, zobs) for yy in yg]),'ro')

        ax13 = fig.add_subplot(3,3,3)
        ax13.plot(zl, self.dHdx(xobs, yobs, zl), 'b-')
        ax13.plot(zl, np.array([ TI.ddx(xobs, yobs, zz) for zz in zl]),'r--')
        ax13.plot(zg, np.array([ TI.ddx(xobs, yobs, zz) for zz in zg]),'ro')

        ax21 = fig.add_subplot(3,3,4)
        ax21.plot(xl, self.dHdy(xl, yobs, zobs), 'b-')
        ax21.plot(xl, np.array([ TI.ddy(xx, yobs, zobs) for xx in xl]),'r--')
        ax21.plot(xg, np.array([ TI.ddy(xx, yobs, zobs) for xx in xg]),'ro')

        ax22 = fig.add_subplot(3,3,5)
        ax22.plot(yl, self.dHdy(xobs, yl, zobs), 'b-')
        ax22.plot(yl, np.array([ TI.ddy(xobs, yy, zobs) for yy in yl]),'r--')
        ax22.plot(yg, np.array([ TI.ddy(xobs, yy, zobs) for yy in yg]),'ro')

        ax23 = fig.add_subplot(3,3,6)
        ax23.plot(zl, self.dHdy(xobs, yobs, zl), 'b-')
        ax23.plot(zl, np.array([ TI.ddy(xobs, yobs, zz) for zz in zl]),'r--')
        ax23.plot(zg, np.array([ TI.ddy(xobs, yobs, zz) for zz in zg]),'ro')

        ax31 = fig.add_subplot(3,3,7)
        ax31.plot(xl, self.dHdtau(xl, yobs, zobs), 'b-')
        ax31.plot(xl, np.array([ TI.ddz(xx, yobs, zobs) for xx in xl]),'r--')
        ax31.plot(xg, np.array([ TI.ddz(xx, yobs, zobs) for xx in xg]),'ro')

        ax32 = fig.add_subplot(3,3,8)
        ax32.plot(yl, self.dHdtau(xobs, yl, zobs), 'b-')
        ax32.plot(yl, np.array([ TI.ddz(xobs, yy, zobs) for yy in yl]),'r--')
        ax32.plot(yg, np.array([ TI.ddz(xobs, yy, zobs) for yy in yg]),'ro')

        ax33 = fig.add_subplot(3,3,9)
        ax33.plot(zl, self.dHdtau(xobs, yobs, zl), 'b-')
        ax33.plot(zl, np.array([ TI.ddz(xobs, yobs, zz) for zz in zl]),'r--')
        ax33.plot(zg, np.array([ TI.ddz(xobs, yobs, zz) for zz in zg]),'ro')
        ax31.set_xlabel('x')
        ax32.set_xlabel('y')
        ax33.set_xlabel('tau')
        ax11.set_ylabel('ddx')
        ax21.set_ylabel('ddy')
        ax31.set_ylabel('ddtau')

#rf = RFMultipole()
#rf.knl = [0, 0, 0, 60, 0, 0*6000]
#rf.frequency = 6*400789598.4106178
#
#rf.init_H()
#rf.report(1.2,1.2,0.5,7,7,15)
#
#plt.show()
