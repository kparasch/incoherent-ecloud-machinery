import numpy as np
from scipy.constants import c as clight
import sympy
import math


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
        print('Hamiltonian is:')
        print(self.Hamiltonian)

    def H(self, x, y, tau):
        return self.Hamiltonian.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdx(self, x, y, tau):
        return self.d_Hamiltonian_dx.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdy(self, x, y, tau):
        return self.d_Hamiltonian_dy.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdtau(self, x, y, tau):
        return self.d_Hamiltonian_dtau.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdxdy(self, x, y, tau):
        return self.d_Hamiltonian_dxdy.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdxdtau(self, x, y, tau):
        return self.d_Hamiltonian_dxdtau.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdydtau(self, x, y, tau):
        return self.d_Hamiltonian_dydtau.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )

    def dHdydtau(self, x, y, tau):
        return self.d_Hamiltonian_dxdydtau.subs(
            [(self.xsymb, x), (self.ysymb, y), (self.tausymb, tau)]
        )


rf = RFMultipole()
rf.knl = [0, 0.5, 1]
rf.frequency = 1

rf.init_H()
