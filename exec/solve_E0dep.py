#! /usr/bin/env python2.7

import vector_potential
import numpy as np
import cmath
import dirac_solver
#import pypar

def heaviside(x, x0=0):

     return 0 if x <= x0 else 1

jxn = []
jyn = []

delta = 0.0
time_end = 40.0

numproc = 1 #pypar.size()
myid = 0 #pypar.rank()
node = 0 #pypar.get_processor_name()

#E0range = np.linspace(40.0, 60.0, 11) * 10**5
E0range = np.array([20.0, 40.0])*10**5
dEproc = len(E0range) / numproc

E0rangep = E0range[myid*dEproc: (myid+1)*dEproc]
if myid == numproc - 1:
   if len(E0rangep)%numproc != 0:
      E0rangep = E0range[myid*dEproc:]

print E0rangep

for E0 in E0rangep:
    print 'E0', E0
    fr = 2.0*np.pi*2.*10**(12)
    A0 = E0/fr
    vf = 10**6
    en = 40*10**5 / fr * 10**6
    en0 = A0 * 10**6
    px0 = en/en0
    #print px0
    prefactor = (1.05*10**(-34)*fr/vf/1.6*10**(19))/ (A0)
    #print prefactor
    #A = vector_potential.CustomEnvelopePulse(
    #      amplitude_E0=1.0, frequency=1.0, envelope = (lambda t: (0 if t>5.7e-13*fr else np.exp(-np.abs(((t-5.7e-13*fr)/2.2e-13/fr))**2.7/1.2)) + (0 if t<=5.7e-13*fr else np.exp(-np.abs(((t-5.7e-13*fr)/1.8e-13/fr))**2.0/3.5))),
    #       Nc=2, cep=1.8)
    A = vector_potential.PulseFromFile(amplitude_E0=1.0, frequency=1.0,
                                       file_name="exec/Adat.dat",
                                       tmax=120.*2.0*np.pi)
    pxlist = np.linspace(-px0, px0, 1)
    pylist = np.linspace(-px0, px0, 1)
    solver = dirac_solver.DiracSolver(A=A, prefactor=prefactor)
    jxn = []
    jyn = []
    for px in pxlist:
        for py in pylist:
            p0 = [px, py]
            phase = np.arccos(p0[0]/np.sqrt(p0[0]**2+p0[1]**2)) #np.arctan(p0[1]/p0[0])
            if p0[1] < 0:
                phase *= -1.0
            #print 'mom', px, py, phase
            pe = -np.sqrt(p0[0]**2 + p0[1]**2 + delta**2)
            norm = np.sqrt((pe-delta)/2.0/(pe))
            phi0 = [norm*np.sqrt(p0[0]**2+p0[1]**2)/(pe-delta)*complex(np.cos(phase/2.0),-np.sin(phase/2.0)),
                    1/np.sqrt(2)*complex(np.cos(phase/2.0), np.sin(phase/2.0))]
            solver.solveTD(t0=time_end-1.0, tf=time_end, dt=0.0001, phi0=phi0, p0=p0, delta=0.0, )
            jxn.append(p0 + solver.jx.tolist())
            jyn.append(p0 + solver.jy.tolist())
    datajx = -np.sum(np.array(jxn), axis=0)[2:]/len(jxn)
    Ei = E0 / 10**5
    np.save('data_total_jx_E0%.2f.npy'%Ei, datajx)
    print 'End'

if myid == 0:
   for i in xrange(1, numproc):
      pypar.receive(source=i)
else:
   pypar.send('', 0)

pypar.finalize()
