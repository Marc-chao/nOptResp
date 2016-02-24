from scipy.integrate import ode
import numpy as np
from math import *
import cmath
from random import randrange, uniform
import utilites

class OpticalSolver:

    def __init__(self, A, prefactor=9.46*10**(-3), tdamp = 10**8):

        self.A = A
        self.prefactor = prefactor
        self.jx = None
        self.jy = None
        self.tarr = None
        self.vec = None
        self.omega_np = None
        self.t_np = None
        self.tdamp = tdamp

    def make_results(self, c, p0):
        jx = []
        jy = []
        n = -np.abs(self.vec[:,0])**2 + np.abs(self.vec[:,1])**2
        rho = self.vec[:,2]
        A_data = np.array([self.A(t)[0] for t in self.tarr])
        omega = utilites.get_Omega(p0=p0,tarr=self.tarr, A=self.A)/self.prefactor
        #for i in xrange(len(c)):
        self.jx = n*(A_data+p0[0])/np.sqrt((p0[0]+A_data)**2+p0[1]**2) \
                                        +(1.j*p0[1]/np.sqrt((p0[0]+A_data)**2+p0[1]**2)*\
                                        (rho*np.exp(-2.j*omega) - np.conjugate(rho)*np.exp(2.j*omega)))
        self.jy =n*(A_data+p0[0])/np.sqrt((p0[0]+A_data)**2+p0[1]**2) \
                                        +(1.j*(p0[0]+A_data)/np.sqrt((p0[0]+A_data)**2+p0[1]**2)*\
                                        (rho*np.exp(-2.j*omega) - np.conjugate(rho)*np.exp(2.j*omega)))
        return None

    def f(self, t, vec0, p0):
        A = self.A
        dt = self.t_np[1] - self.t_np[0]
        index = np.where(self.t_np <= t)[0][-1]
        try:
            omega = self.omega_np[index] + np.sqrt((p0[0]+A(t)[0])**2 +
                                                                            (p0[1]+A(t)[1])**2)*(t-self.t_np[index])
        except:
            print t
            omega = self.omega_np[-1]
        tetta_0 = np.sign(p0[1])*\
                       np.arccos((p0[0]+A(t)[0])/np.sqrt((p0[0]+A(t)[0])**2 +
                                                                            (p0[1]+A(t)[1])**2))
        tetta_1 = np.sign(p0[1])*\
                        np.arccos((p0[0]+A(t+dt)[0])/np.sqrt((p0[0]+A(t+dt)[0])**2+
                                                                                   (p0[1]+A(t+dt)[1])**2))
        tetta_dot = (tetta_1 - tetta_0)/dt
        a = -0.5j * tetta_dot * (np.conjugate(vec0[2]) * np.exp(2.0j*omega/self.prefactor) -\
                vec0[2] * np.exp(-2.0j*omega/self.prefactor))
        b = 0.5j * tetta_dot * (np.conjugate(vec0[2]) * np.exp(2.0j*omega/self.prefactor) -\
                vec0[2] * np.exp(-2.0j*omega/self.prefactor))
        c = -0.5j * tetta_dot * (vec0[1]-vec0[0]) * np.exp(2.0j*omega/self.prefactor) - vec0[2]/self.tdamp
        return np.array([a,  b, c])

    def solveTD(self, t0=0.0, tf=1.0, dt=0.001,
                       vec0=[complex(1.0,0.0), complex(0.0,0.0), complex(0.0,0.0)],
                       p0=[-0.75, 0.02]):


        self.t_np = np.arange(0, 1.2*tf, dt/10.0)
        self.omega_np = utilites.get_Omega(p0=p0, tarr = self.t_np, A=self.A)

        r = ode(self.f).set_integrator('zvode', method='bdf', order=2)
        r.set_initial_value(vec0, t0).set_f_params(p0).set_jac_params(p0)
        vec = []
        tarr = []
        vec.append(vec0)
        tarr.append(0.0)

        while r.successful() and r.t < tf:
            r.integrate(r.t+dt)
            vec.append(r.y)
            tarr.append(r.t)

        self.tarr = np.array(tarr)
        self.vec = np.array(vec)
        self.make_results(c=vec, p0=p0)
        #return phi
