import unittest
import vector_potential
import numpy as np
import cmath
from scipy.ndimage.filters import gaussian_filter
import optical_equation_solver
import utilites

class TestEquations(unittest.TestCase):

    def setUp(self):
        pass

    def test_optical_equations(self):
        E0 = 40.0*10**(5)
        fr = 2.0*np.pi*2.*10**(12)
        A0 = E0/fr
        vf = 10**6
        prefactor = (1.05*10**(-34)*fr/vf/1.6*10**(19))/ (E0/fr)
        A = vector_potential.GaussianEnvelopePulse(amplitude_E0=1.0, frequency=1.0,
                                           Nc=2, cep=0.0, tc=10.0)

        print "PASS"
