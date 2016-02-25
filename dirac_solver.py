from scipy.integrate import ode
import numpy as np
from math import *

class DiracSolver:

    def __init__(self, A, prefactor=9.46*10**(-3)):

        self.A = A
        self.prefactor = prefactor
        self.jx = None
        self.jy = None
        self.tarr = None
        self.phi = None

    def make_results(self, c):
        #jx = []
        #jy = []

        #for i in xrange(len(c)):
         #   jx.append(np.conjugate(c[i][0]) * c[i][1] + np.conjugate(c[i][1]) * c[i][0])
          #  jy.append(1j*(-np.conjugate(c[i][0]) * c[i][1] + np.conjugate(c[i][1]) * c[i][0]))
        self.jx = np.conjugate(c[:, 0]) * c[:, 1] + np.conjugate(c[:, 1]) * c[:, 0] #np.array(jx)
        self.jy = 1.j*(-np.conjugate(c[:, 0]) * c[:,1] + np.conjugate(c[:, 1]) * c[:, 0]) #np.array(jy)
        return None

    def f(self, t, phi0, p0, delta, ksi, ksi0):

        A = self.A

        a = 1/(1j * self.prefactor) * (delta * phi0[0] + \
                                                     (p0[0]+A(t)[0]-1.0j*(p0[1]+A(t)[1]))*phi0[1] -\
            ksi * ((p0[0]+A(t)[0])**2 - (p0[1]+A(t)[1])**2 +
                      2.0j*(p0[0]+A(t)[0])*(p0[1]+A(t)[1]))* phi0[1]) -\
            ksi0 * ((p0[0]+A(t)[0])**2 + (p0[1]+A(t)[1])**2) * phi0[0]
        b = 1/(1j * self.prefactor) * (-delta * phi0[1] + \
                                                     (p0[0]+A(t)[0]+1.0j*(p0[1]+A(t)[1]))*phi0[0] -\
            ksi * ((p0[0]+A(t)[0])**2 - (p0[1]+A(t)[1])**2 -
                      2.0j*(p0[0]+A(t)[0])*(p0[1]+A(t)[1]))*phi0[0]) -\
            ksi0 * ((p0[0]+A(t)[0])**2 + (p0[1]+A(t)[1])**2) * phi0[1]

        return np.array([a,  b])

    def solveTD(self, t0=0.0, tf=1.0, dt=0.0001,
                       phi0=[complex(0.0,0.0), complex(1.0,0.0)],
                       p0=[-0.75, 0.02], delta=0.0, ksi=0.0, ksi0=0.0):

        self.phase_old = 0.0
        self.told = 0.0
        r = ode(self.f).set_integrator('zvode', method='bdf', order=2)
        r.set_initial_value(phi0, t0).set_f_params(p0, delta, ksi, ksi0).set_jac_params(p0)
        phi = []
        tarr = []
        phi.append(phi0)
        tarr.append(0.0)
        i=0
        while r.successful() and r.t < tf:
            i += 1
            r.integrate(r.t+dt)
            if np.mod(i, 1.0) == 0.0:
                phi.append(r.y)
                tarr.append(r.t)
        self.tarr = np.array(tarr)
        self.phi = np.array(phi)
        self.make_results(c=self.phi)
        #return phi
