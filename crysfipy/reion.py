import crysfipy.const as C
from crysfipy.const import ion
from crysfipy.cfmatrix import *
import numpy as np
from numpy import diag, conj, transpose, dot
from numpy.linalg import eig
import numbers

class cfpars:
    """Set of crystal field parameters"""
    
    pars = {
    "c": ["B40", "B60"],
    "h": ["B20", "B40", "B44", "B66"],
    "t": ["B20", "B40", "B44", "B60", "B64"],
    "o": ["B20", "B22", "B40", "B42", "B44", "B60", "B62", "B64", "B66"],
    }
    
    def __init__(self, *args, **kwargs):
        skipNonnamed = False
        self.sym = ""  #orthorombic or more general symetry (=no constrains)
        #clear all params
        for name in self.pars["o"]: self.asignParameter(name, 0)
        #check named args
        for name, value in kwargs.items():
            #check symetry constrain
            if name == "sym":
                self.sym = self.asignSymmetry(value)    
            if name[0] == "B":  #it is a crystal field
                skipNonnamed = True
                self.asignParameter(name, value)
        #check non named args
        if self.sym == "":
            for value in args:
                if isinstance(value, str):
                    self.sym = self.asignSymmetry(value)
                    break  #just consider first string in list
        if self.sym == "": self.sym = "o" #default
        if not skipNonnamed:
            #read
            i = 0
            for value in args:
                if isinstance(value, numbers.Real):
                    if i >= len(self.pars[self.sym]): break
                    self.asignParameter(self.pars[self.sym][i], value)
                    i+=1
        #do symmetry magic
        if self.sym == "c":
            raise NotImplementedError("Cubic symmetry not yet supported!")
            
    def __str__(self):
        ret = ""
        for name in self.pars["o"]: ret += self.printParameter(name) + "\n"
        return ret
    
            
    def asignSymmetry(self, inp):
        if inp[0] == "t" or  inp[0] == "h" or  inp[0] == "c":
            return inp[0]
        return ""
                  
    def asignParameter(self, name, value):
        if name == "B20": self.B20 = value
        if name == "B22": self.B22 = value
        if name == "B40": self.B40 = value
        if name == "B42": self.B42 = value
        if name == "B44": self.B44 = value
        if name == "B60": self.B60 = value
        if name == "B62": self.B62 = value
        if name == "B64": self.B64 = value
        if name == "B66": self.B66 = value
        
    def printParameter(self, name):
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
    """Object representing rare-earth ion in CF potential"""

    def __init__(self, name, field, cfp):
        self.name = name
        self.field = field
        if type(cfp) is list:
            cfp = cfpars(*cfp)
        self.cfp = cfp
        
        self.H = np.array(field)
        self.H_size = np.sqrt(dot(self.H, self.H.conj().transpose()))
        if self.H_size > 0:
            self.H_direction = self.H / self.H_size
        self.calculate()
    
    
    
    def getlevels(self):
        """Calculate degeneracy of the levels and sort the matrix"""
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
      
    def calculate(self):
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
        
        
        hamiltonian =  \
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
        
        E, U = eig(hamiltonian);
        
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



def rawneutronint(E, deg, J2, gJ, T):
    """Returns transition intensities in barn."""
    """E - matrix of energy levels in meV"""
    """J2 - matrix of squared J"""
    r02 = const.R0 * C.R0  *1e28 # to have value in barn
    c = np.pi * r02 * gJ * gJ
    
    prst = np.exp(-E*C.eV2K/T)
    Z = sum(prst*deg) #multiply with degeneracy of the level
    prst = prst / Z
    trans_int = J2 * prst[:, np.newaxis] * c  #transition intensities in barn
    return trans_int

def neutronint(ion, T, direction = "t"):
    """Returns matrix of energy and transition intensity at given temperature"""
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
    tint = rawneutronint(ion.deg_e[:,0], ion.deg_e[:,1], J2, ion.gJ, T).flatten()
    tint = tint[jumps.argsort()]
    return np.array([np.sort(jumps),tint]) 

def rawsusceptibility(energy, moment, H_direction, H_size, T):
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
            y.append(rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, temp))
    except TypeError as te:
        y = rawsusceptibility(ion.energy, ion.moment, ion.H_direction, ion.H_size, T)
    return y