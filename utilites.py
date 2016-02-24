import numpy as np
import matplotlib.pylab as plt
import math
try:
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
except:
    print 'Warning(vector_potentian): no module matplotlib'
    pass

#Does not work
#def energy(p0, A, tarr, delta=0, ksi=0, vf = 10**6,):

#    return vf * A0 * (np.sqrt(parr[0]**2+parr[1]**2)) -\
 #            ksi * A0**2 * (parr[0]**2+parr[1]**2) * \
  #           np.cos(3.0*np.arccos(parr[0]/np.sqrt(parr[1]**2 + parr[0]**2)))

def energy(parr, vF, A0, ksi, ksi0 = 0, delta=0):
    return np.sqrt(delta**2 + A0**2 * np.sqrt(parr[0]**2+parr[1]**2)**2 *\
                (vF**2 +\
                ksi**2 * A0**2 * (parr[0]**2+parr[1]**2) -\
                2.0 *  ksi * vF * A0 * np.sqrt(parr[0]**2+parr[1]**2) *\
                    np.cos(3.0*np.arccos(parr[0]/np.sqrt(parr[1]**2 + parr[0]**2)))))


def intraband_current(p0, A, tarr, delta=0.0):
    return np.array([-(p0[0]+A(t)[0])/np.sqrt((p0[0]+A(t)[0])**2 +
                              (p0[1]+A(t)[1])**2 + delta**2) for t in tarr])

def intraband_current_y(p0, A, tarr, delta=0.0):

    return np.array([(p0[1]+A(t)[1])/np.sqrt((p0[0]+A(t)[0])**2 +
                              (p0[1]+A(t)[1])**2 + delta**2) for t in tarr])

# In case of warping these inerband current is only an approximation.
# TO DO: make a correct formula
def intraband_current_war(p0, A, tarr, ksi0=0.0, ksi=0.0):

    return np.array([-(p0[0] + A(t)[0])/np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)\
            + 2.0 * np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)*(p0[0] + A(t)[0])*(ksi0 +
                    ksi*np.cos(3.0*np.arccos((p0[0] + A(t)[0])/np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2))))\
            + 3.0 * ksi * np.abs((p0[1] + A(t)[1])) * \
                    np.sin(3.0*np.arccos((p0[0] + A(t)[0])  / np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)))\
            for t in tarr])

def intraband_current_war_y(p0, A, tarr, ksi0=0.0, ksi=0.0):
   return ([-(p0[1] + A(t)[1])/np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)\
            + 2.0 * np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)*(p0[0] + A(t)[0])*(ksi0 +
                    ksi*np.cos(3.0*np.arccos((p0[0] + A(t)[0])/np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2))))\
            + 3.0 * ksi * (p0[0] + A(t)[0])*np.sign(p0[1] + A(t)[1]) * \
                    np.sin(3.0*np.arccos((p0[0] + A(t)[0])  / np.sqrt((p0[1] + A(t)[1])**2 + (p0[0] + A(t)[0])**2)))\
            for t in tarr])


def state_down(p0, A, tarr, prefactor, omega=None, delta=0):

    if omega is None:
       omega = get_Omega(p0=p0, A=A, tarr=tarr)

    A_0 = np.array([A(t)[0] for t in tarr])
    A_1 = np.array([A(t)[1] for t in tarr])

    phase = np.arccos((p0[0]+A_0)/np.sqrt((p0[0]+A_0)**2+(p0[1]+A_1)**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0

    pe = -np.sqrt((p0[0]+A_0)**2 + (p0[1]+A_1)**2 + delta**2)
    norm = np.sqrt((pe-delta)/2.0/(pe))
    phi0m = np.array([np.sqrt((p0[0]+A_0)**2+(p0[1]+A_1)**2)/(pe-delta)*(np.cos(phase/2.0)-1.j*np.sin(phase/2.0)),
    np.cos(phase/2.0) + 1.j*np.sin(phase/2.0)]) * norm * (np.cos(omega/prefactor) + 1.j*np.sin(omega/prefactor))
    return phi0m.T

def state_up(p0, A, tarr, prefactor, omega=None, delta=0):

    if omega is None:
       omega = get_Omega(p0=p0, A=A, tarr=tarr)

    A_0 = np.array([A(t)[0] for t in tarr])
    A_1 = np.array([A(t)[1] for t in tarr])

    phase = np.arccos((p0[0]+A_0)/np.sqrt((p0[0]+A_0)**2+(p0[1]+A_1)**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0

    pe = np.sqrt((p0[0]+A_0)**2 + (p0[1]+A_1)**2 + delta**2)
    norm = np.sqrt((pe-delta)/2.0/(pe))
    phi0u = np.array([np.sqrt((p0[0]+A_0)**2+(p0[1]+A_1)**2)/(pe-delta)*(np.cos(phase/2.0)-1.j*np.sin(phase/2.0)),
    (np.cos(phase/2.0)+1.j*np.sin(phase/2.0))])*norm*(np.cos(omega/prefactor)-1.j*np.sin(omega/prefactor))
    return phi0u.T


def find_Coeffs(state, state_up, state_down):

    C_down = [np.dot(np.conj(state_down[i]).T, state[i]) for i in xrange(len(state))]
    C_up = [np.dot(np.conj(state_up[i]).T, state[i]) for i in xrange(len(state))]
    return C_down, C_up

# Initial wave functions
def state_down0(p0, delta=0):
    phase = np.arccos((p0[0])/np.sqrt((p0[0])**2+(p0[1])**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0 
    pe = -np.sqrt((p0[0])**2 + (p0[1])**2 + delta**2)
    norm = np.sqrt((pe-delta)/2.0/pe)
    phi0m = norm*np.array([np.sqrt((p0[0])**2+(p0[1])**2)/(pe-delta)*np.exp(-0.5j*phase),
                                           np.exp(0.5j*phase)])
    return phi0m


def state_up0(p0, delta=0):
    phase = np.arccos((p0[0])/np.sqrt((p0[0])**2+(p0[1])**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0
    pe = np.sqrt((p0[0])**2 + (p0[1])**2 + delta**2)
    norm = np.sqrt((pe-delta)/2.0/pe)
    phi0u = norm*np.array([np.sqrt((p0[0])**2+(p0[1])**2)/(pe-delta)*np.exp(-0.5j*phase),
                                           np.exp(0.5j*phase)])
    return phi0u


# Wave function in the case of warping
# TO DO: add a gap.
def state_down0_war(p0, vF, ksi=0.0):
    phase = np.arccos((p0[0])/np.sqrt((p0[0])**2+(p0[1])**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0 
    pe = np.sqrt(p0[0]**2 + p0[1]**2)

    phi0m = 1.0/np.sqrt(2.0)*\
            np.array([np.sqrt((vF*np.exp(-1.j*phase) - ksi * pe * np.exp(2.j*phase))/ \
                               (vF*np.exp(1.j*phase) - ksi * pe * np.exp(-2.j*phase))),1.0])
    return phi0m

def state_up0_war(p0, vF, ksi=0.0):
    phase = np.arccos((p0[0])/np.sqrt((p0[0])**2+(p0[1])**2))  #np.arctan(p0[1]/p0[0])
    if p0[1] < 0:
        phase *= -1.0
    pe = np.sqrt(p0[0]**2 + p0[1]**2)

    phi0u = 1.0/np.sqrt(2.0)*\
            np.array([-np.sqrt((vF*np.exp(-1.j*phase) - ksi * pe * np.exp(2.j*phase))/ \
                               (vF*np.exp(1.j*phase) - ksi * pe * np.exp(-2.j*phase))),1.0])
    return phi0u

# Have to be corrected: WRONG!
#def state_down0_war_gap(p0, vF=vf, ksi=0.0, delta=0.0):
#    phase = np.arccos((p0[0])/np.sqrt((p0[0])**2+(p0[1])**2))  #np.arctan(p0[1]/p0[0])
#   if p0[1] < 0:
#        phase *= -1.0 
#    pe = np.sqrt(p0[0]**2 + p0[1]**2)
#    norm = (delta - np.sqrt(delta**2 + pe**2 * (vF*np.exp(-1.j*phase) - ksi * pe * np.exp(2.j*phase))*\
#                                      (vF*np.exp(1.j*phase) - ksi * pe * np.exp(-2.j*phase))))/\
#            np.sqrt((vF*pe*np.exp(-1.j*phase) - ksi * pe**2 * np.exp(2.j*phase))*\
#                    (vF*pe*np.exp(1.j*phase) - ksi * pe**2 * np.exp(-2.j*phase)))
#    phi0m = norm *\
#            np.array([(vF*pe*np.exp(-1.j*phase) - ksi * pe**2 * np.exp(2.j*phase))/ \
#                               (delta - np.sqrt(delta**2 + pe**2 * (vF*np.exp(-1.j*phase) - ksi * pe * np.exp(2.j*phase))*\
#                                      (vF*np.exp(1.j*phase) - ksi * pe * np.exp(-2.j*phase)))),1.0])
#    return phi0m


def get_Omega(p0, tarr, A):
    om=[]
    om0 = 0.0
    om.append(om0)
    dt = tarr[-1] - tarr[-2]
    for j in xrange(1, len(tarr)):
        om0 += np.sqrt((p0[0]+A(tarr[j])[0])**2 + (p0[1]+A(tarr[j])[1])**2)*dt
        om.append(om0)
    return np.array(om)

def window_function(t, width=0.10):
    if np.abs(t)<=1.0/width:
        return (1.0 + np.cos(width*np.pi*t))/2.0
    else: return 0.0

def windowed_fourier_transform(func, twft=np.linspace(0.0,20.0,100.0),
                                                  tarr=np.linspace(0.0,20.0,100.0),
                                                  width=1.5):
    wft = []
    for ts in twft:
        wf = np.array([window_function(t=t-ts, width=width) for t in tarr])
        a = np.fft.fft(func*wf)
        wft.append(a**2)
    return wft

