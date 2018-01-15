# Copyright 2014-2018 Petr Čermák, Jan Zubáč and Karel Pajskr
# This file is part of CrysFiPy.
# CrysFiPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# CrysFiPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# <http://www.gnu.org/licenses/>.


import crysfipy.const as C
from crysfipy.const import ion
from crysfipy.cfmatrix import *
import numpy as np
from numpy import diag, conj, transpose, dot
from numpy.linalg import eig
import numbers

class cfpars:
    """Class representing set of crystal field parameters.

    It simplifies the creation of the CF parameter sets considering symmetry
    of the environment.
    Other modules expect that CF parameters are in *meV* (SI units). 
    But if you just want to diagonalize Hamiltonian, it is possible to use *K* (and results will be in *K*).

    Initialization can be done with named arguments or without. 
    If arguments are not named, symmetry is considered from the first string argument.
    Stevens parameters are considered differently for different symmetries in following order:

    | cubic: B40, B60
    | hexagonal: B20, B40, B44, B66
    | tetragonal: B20, B40, B44, B60, B64
    | orthorombic: B20, B22, B40, B42, B44, B60, B62, B64, B66
    

    Attributes:
        BXY (float, optional): Attribute corresponding to :math:`B_X^Y` Stevens Parameters.
            See Hutchings_.  If at least one CF parameter is specified as a named argument, 
            non-named numerical parameters are ignored.
        sym (str, optional): Symmetry of the crystal field

            | c - cubix
            | h - hexagonal
            | t - tetragonal
            | o - orthorombic (default)
           
    Examples:
        Create set of CF parameters by named parameters:

        >>> print(cfpars(sym = "c", B40 = 10))
        Set of CF parameters for cubic symmetry:
        B40 = 10.0000
        B60 = 0.0000
        B44 = 50.0000
        B64 = 0.0000

        Use of non-named parameters:

        >>> print(cfpars("c", 10, 1))
        Set of CF parameters for cubic symmetry:
        B40 = 10.0000
        B60 = 1.0000
        B44 = 50.0000
        B64 = -21.0000

    """
    
    pars = {
    "c": ["cubic", ["B40", "B60", "B44", "B64"]],            
    "h": ["hexagonal", ["B20", "B40", "B44", "B66"]],              
    "t": ["tetragonal", ["B20", "B40", "B44", "B60", "B64"]],            
    "o": ["orthorombic", ["B20", "B22", "B40", "B42", "B44", "B60", "B62", "B64", "B66"]],  
    }
    
    def __init__(self, *args, **kwargs):
        skipNonnamed = False
        self.sym = ""  #orthorombic symetry (=no constrains)
        #clear all params
        for name in self.pars["o"][1]: self._asignParameter(name, 0)
        #check named args
        for name, value in kwargs.items():
            #check symetry constrain
            if name == "sym":
                self.sym = self._asignSymmetry(value)    
            if name[0] == "B":  #it is a crystal field
                skipNonnamed = True
                self._asignParameter(name, value)
        #check non named args
        if self.sym == "":
            for value in args:
                if isinstance(value, str):
                    self.sym = self._asignSymmetry(value)
                    break  #just consider first string in list
        if self.sym == "": self.sym = "o" #default
        if not skipNonnamed:
            #read
            i = 0
            for value in args:
                if isinstance(value, numbers.Real):
                    if i >= len(self.pars[self.sym][1]): break
                    self._asignParameter(self.pars[self.sym][1][i], value)
                    i+=1
        #do symmetry magic
        if self.sym == "c":
            self.B44 = 5 * self.B40;
            self.B64 = -21 * self.B60;
            
    def __str__(self):
        ret = "Set of CF parameters for %s symmetry:" % (self.pars[self.sym][0])
        for name in self.pars[self.sym][1]: ret += self._printParameter(name) + "\n"
        return ret
    
            
    def _asignSymmetry(self, inp):
        if inp[0] == "t" or  inp[0] == "h" or  inp[0] == "c":
            return inp[0]
        return ""
                  
    def _asignParameter(self, name, value):
        if name == "B20": self.B20 = value
        if name == "B22": self.B22 = value
        if name == "B40": self.B40 = value
        if name == "B42": self.B42 = value
        if name == "B44": self.B44 = value
        if name == "B60": self.B60 = value
        if name == "B62": self.B62 = value
        if name == "B64": self.B64 = value
        if name == "B66": self.B66 = value
        
    def _printParameter(self, name):
        if name == "B20": return "%s = %.4f" % (name, self.B20)
        if name == "B22": return "%s = %.4f" % (name, self.B22)
        if name == "B40": return "%s = %.4f" % (name, self.B40)
        if name == "B42": return "%s = %.4f" % (name, self.B42)
        if name == "B44": return "%s = %.4f" % (name, self.B44)
        if name == "B60": return "%s = %.4f" % (name, self.B60)
        if name == "B62": return "%s = %.4f" % (name, self.B62)
        if name == "B64": return "%s = %.4f" % (name, self.B64)
        if name == "B66": return "%s = %.4f" % (name, self.B66)




class re:
    """Object representing rare-earth ion in CF potential
    
    Attributes:
        name (str): Name of the ion.
        field (1D array of floats): external magnetic field applied in *T*.
        cfp (:obj:`crysfipy.reion.cfpars`): Crystal field parameters
        calculate (bool, optional): If true (default) then it automatically diagonalizes
            Hamiltonian and calculates energy levels.

    Examples:
        
        >>> ce = re("Ce", [0,0,0], ["c", 10])
        >>> print(ce)
        Energy levels:
        E(0) =	0.0000	 2fold-degenerated
        E(1) =	3600.0000	 4fold-degenerated
    """

    def __init__(self, name, field, cfp, calculate = True):
        self.name = name
        self.field = field
        if type(cfp) is list:
            cfp = cfpars(*cfp)
        self.cfp = cfp
        
        self.H = np.array(field)
        self.H_size = np.sqrt(dot(self.H, self.H.conj().transpose()))
        if self.H_size > 0:
            self.H_direction = self.H / self.H_size
        if (calculate):
            self.getlevels()
    
    
    
    def getlevels(self):
        """Calculate degeneracy of the levels and sort the matrix"""

        self._calculate()
        #orthogonalization not needed?
        #self.rawev,R = np.linalg.qr(self.rawev)
        
        #change the sign to be positive :)
        self.rawev = self.rawev * np.sign(np.sum(self.rawev, axis=0))

        self.energy = self.rawenergy - min(self.rawenergy)     # shift to zero level	
        #sorting
        self.ev = self.rawev[:,self.rawenergy.argsort()]
        self.energy = np.sort(self.energy)
        #get sorted Jx,Jz and Jz 
        self.Jz = dot(dot(self.ev.conj().transpose(), J_z(self.J)), self.ev)
        self.Jx = dot(dot(self.ev.conj().transpose(), J_x(self.J)), self.ev)
        self.Jy = dot(dot(self.ev.conj().transpose(), J_y(self.J)), self.ev)
        #calculate J^2 matrices
        self.Jx2 = np.square(np.abs(self.Jx))
        self.Jy2 = np.square(np.abs(self.Jy))
        self.Jz2 = np.square(np.abs(self.Jz))

        #calculate degeneracy
        deg_e = []
        levels = 0
        deg_e.append([self.energy[0], 0])
        for x in self.energy:
            if not np.isclose(deg_e[levels][0], x):
                levels+=1
                deg_e.append([x, 1])
            else:
                deg_e[levels][1] += 1
        
        levels += 1  # started at zero
        
        #empty degenerate level transition matrices
        self.deg_Jx2 = np.zeros((levels,levels))
        self.deg_Jy2 = np.zeros((levels,levels))
        self.deg_Jz2 = np.zeros((levels,levels))
        #sum degenerated levels
        u = 0
        v = 0
        for i, x in enumerate(deg_e):
            v = 0
            for j, y in enumerate(deg_e):
                self.deg_Jx2[i,j] = np.sum(self.Jx2[u:u+x[1],v:v+y[1]])
                self.deg_Jy2[i,j] = np.sum(self.Jy2[u:u+x[1],v:v+y[1]])
                self.deg_Jz2[i,j] = np.sum(self.Jz2[u:u+x[1],v:v+y[1]])
                v+=y[1]
            u+=x[1]        
            
        #calculate Jt2 for polycrystal
        self.deg_Jt2 = 2.0 / 3 * (self.deg_Jx2 + self.deg_Jy2 + self.deg_Jz2)
        self.deg_e = np.array(deg_e)   
        
    def __str__(self):
        """Nice printout of calculated parameters"""
        ret = ""
        ret += "Energy levels:\n"
        for i,x in enumerate(self.deg_e):
            ret += "E({:d}) =\t{:.4f}\t{:2d}fold-degenerated\n".format(i, x[0], int(x[1]))
        return ret
      
    def _calculate(self):
        """Calculates energy splitting in CF potential"""

        i = ion(self.name)
        self.J = i.J
        self.p1 = np.ones((i.J2p1,1), float);      # column vector of ones
        self.Jx = J_x(i.J); 
        self.Jy = J_y(i.J);
        self.Jz = J_z(i.J);                        # prepare matrices
        self.gJ = i.gJ
        self.moment = np.zeros((i.J2p1,3), float)  # matrix with projection of moments to x, y, z directions for all J2p1 levels
        
        B = self.cfp
        
        
        self.hamiltonian =  \
            B.B20 * O_20(i.J) + \
            B.B22 * O_22(i.J) + \
            B.B40 * O_40(i.J) + \
            B.B42 * O_42(i.J) + \
            B.B44 * O_44(i.J) + \
            B.B60 * O_60(i.J) + \
            B.B62 * O_62(i.J) + \
            B.B64 * O_64(i.J) + \
            B.B66 * O_66(i.J) + \
            C.uB * i.gJ * (self.Jx * self.H[0] + self.Jy * self.H[1] + self.Jz * self.H[2])
        
        E, U = eig(self.hamiltonian);
        
        self.rawenergy = np.real(E);  
        if sum(np.iscomplex(E)) > 0:
            raise ValueError('Final energies are complex!')
        
        self.Jz = dot(dot(U.conj().transpose(), self.Jz), U)   # conversion of matrices to the basis of eigenvectors
        self.Jx = dot(dot(U.conj().transpose(), self.Jx), U)   # it is then easier to calculate <i|J|j>
        self.Jy = dot(dot(U.conj().transpose(), self.Jy), U)
        
        self.moment[:,0] = - i.gJ * np.real(diag(self.Jx))     # calculation of projections of moments for every eigenvector   
        self.moment[:,1] = - i.gJ * np.real(diag(self.Jy))                      
        self.moment[:,2] = - i.gJ * np.real(diag(self.Jz))
        
        self.rawev = U



def _rawneutronint(E, deg, J2, gJ, T):
    """Returns transition intensities in barn.

    Args:
        E (1D array of floats): matrix of energy levels in meV
        deg (1D array of floats): degeneracies of the levels
        J2 (2D array of floats): matrix of squared J
        gJ (float): Landé factor
        T (float): temperature in *K*
    """
    r02 = C.R0 * C.R0  *1e28 # to have value in barn
    c = np.pi * r02 * gJ * gJ
    
    prst = np.exp(-E*C.eV2K/T)
    Z = sum(prst*deg) #multiply with degeneracy of the level
    prst = prst / Z
    trans_int = J2 * prst[:, np.newaxis] * c  #transition intensities in barn
    return trans_int

def neutronint(ion, T, direction = "t"):
    """Returns matrix of energy and transition intensity at given temperature
    
    Args:
        ion (:obj:`crysfipy.reion.re`): Rare-earth ion object
        T (float): temperature in *K*
        direction (str): Direction of the Q in which to calculate
            
            | t - powder (default)
            | x - using :math:`J_x`
            | y - using :math:`J_y`
            | z - using :math:`J_z`
        
    """
    jumps = ion.deg_e[:,0] - ion.deg_e[:,0][:, np.newaxis]
    jumps = jumps.flatten()
    if direction == "x":
        J2 = ion.deg_Jx2
    elif direction == "y":
        J2 = ion.deg_Jy2
    elif direction == "z":
        J2 = ion.deg_Jz2
    else:
        J2 = ion.deg_Jt2
    tint = _rawneutronint(ion.deg_e[:,0], ion.deg_e[:,1], J2, ion.gJ, T).flatten()
    tint = tint[jumps.argsort()]
    return np.array([np.sort(jumps),tint]) 

def _rawsusceptibility(energy, moment, H_direction, H_size, T):
    """Returns susceptibility calculated for energy levels at given temperature"""

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
            y.append(_rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, temp))
    except TypeError as te:
        y = _rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, T)
    return y