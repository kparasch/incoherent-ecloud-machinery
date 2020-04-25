import numpy as np
from scipy.constants import c, e, m_p 

from contextlib import contextmanager
import sys, os

from PyHEADTAIL.machines.synchrotron import Synchrotron

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def sigmas(epsn_x=1.7e-6,epsn_y=1.7e-6, sigma_z=0.09, 
                  beta_x=None, beta_y=None, D_x=None, D_y=None,
                  alpha_mom_compaction=None, circumference=None,
                  rf_harmonic=None, V_rf=None, gamma0=None):

    beta0 = np.sqrt(1-1./gamma0**2)
    n_macroparticles = 3000000

    with suppress_stdout():
        machine = Synchrotron(optics_mode='smooth', charge=e,
                    mass=m_p, p0=beta0*gamma0*m_p*c, n_segments=1, circumference=circumference,
                    beta_x=beta_x, D_x=D_x, beta_y=beta_y, D_y=D_y, accQ_x=0.27,
                    accQ_y=0.29, Qp_x=20., Qp_y=20., longitudinal_mode='non-linear',# Q_s=Q_s,
                    alpha_mom_compaction=alpha_mom_compaction, h_RF=rf_harmonic, V_RF=V_rf,
                    dphi_RF=0., p_increment=0., RF_at='end_of_transverse')

        Q_s = machine.Q_s

        machine_linear = Synchrotron(optics_mode='smooth', charge=e,
                    mass=m_p, p0=beta0*gamma0*m_p*c, n_segments=1, circumference=circumference,
                    beta_x=beta_x, D_x=D_x, beta_y=beta_y, D_y=D_y, accQ_x=0.27,
                    accQ_y=0.29, Qp_x=20., Qp_y=20., longitudinal_mode='linear', Q_s=Q_s,
                    alpha_mom_compaction=alpha_mom_compaction, h_RF=rf_harmonic,# V_RF=V_rf,
                    dphi_RF=0., p_increment=0., RF_at='end_of_transverse')

        bunch_nonlinear = machine.generate_6D_Gaussian_bunch_matched(n_macroparticles=n_macroparticles, intensity=1,
                                                 epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)
        bunch_linear = machine_linear.generate_6D_Gaussian_bunch(n_macroparticles=n_macroparticles, intensity=1,
                                                 epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)

    sigma_x_nonlinear = bunch_nonlinear.sigma_x()
    sigma_y_nonlinear = bunch_nonlinear.sigma_y()
    sigma_x_linear = bunch_linear.sigma_x()
    sigma_y_linear = bunch_linear.sigma_y()

    return sigma_x_nonlinear, sigma_y_nonlinear, sigma_x_linear, sigma_y_linear

