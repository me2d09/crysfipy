import const
from const import ion
from cfmatrix import *
import numpy as np
from numpy import diag, conj, transpose, dot
from numpy.linalg import eig


class re:
    """Object representing rare-earth ion in CF potential"""

    def __init__(self, name, field, cfpars):
        self.name = name
        self.field = field
        self.cfpars = cfpars
        
        self.H = np.array(field)
        self.H_size = np.sqrt(dot(self.H, self.H.conj().transpose()))
        self.H_direction = self.H / self.H_size
        self.calculate()
        
    def calculate(self):
        """Calculates energy splitting in CF potential"""

        i = ion(self.name)
        self.p1 = np.ones((i.J2p1,1), float);      # column vector of ones
        self.Jx = J_x(i.J); 
        self.Jy = J_y(i.J);
        self.Jz = J_z(i.J);                        # prepare matrices
        self.moment = np.zeros((i.J2p1,3), float)  # matrix with projection of moments to x, y, z directions for all J2p1 levels
        
        B = self.cfpars
        
        
        hamiltonian = B[0] * O_20(i.J) + B[1] * O_40(i.J) + B[2] * O_44(i.J) + B[3] * O_60(i.J) + B[4] * O_64(i.J) + \
            const.uB * i.gJ * (self.Jx * self.H[0] + self.Jy * self.H[1] + self.Jz * self.H[2])
        
        E, U = eig(hamiltonian);
        
        
        if sum(np.iscomplex(E)) > 0:
            raise ValueError('Final energies are complex!')
        self.rawenergy = np.real(E);                           
        self.energy = self.rawenergy - min(self.rawenergy)     # shift to zero level	
        
        self.Jz = dot(dot(U.conj().transpose(), self.Jz), U)   # conversion of matrices to the basis of eigenvectors
        self.Jx = dot(dot(U.conj().transpose(), self.Jx), U)   # it is then easier to calculate <i|J|j>
        self.Jy = dot(dot(U.conj().transpose(), self.Jy), U)
        
        self.moment[:,0] = - i.gJ * np.real(diag(self.Jx))     # calculation of projections of moments for every eigenvector   
        self.moment[:,1] = - i.gJ * np.real(diag(self.Jy))                      
        self.moment[:,2] = - i.gJ * np.real(diag(self.Jz))


def rawsusceptibility(energy, moment, H_direction, H_size, T):
    """Returns susceptibility calculated for ënergy levels at given temperature"""

    prst = np.exp(-energy/T)
    Z = sum(prst);                                    # canonical partition function 
    prst = prst / Z;
    overal_moment = dot(prst, moment);
    return dot(overal_moment, H_direction.conj().transpose()) / H_size

def susceptibility(ion, T):
    """Returns susceptibility calculated for given ion at given temperature"""

    y = []
    try:
        for temp in T:
            y.append(rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, temp))
    except TypeError, te:
        y = rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, T)
    return y