from xlpython import *
from CrysFiPy.reion import re, rawsusceptibility as susc
import CrysFiPy.const as C
import numpy as np

@xlfunc
def GetJ2P1(name):
    '''Get number of CF energy levels'''
    return C.ion(name).J2p1


@xlfunc
def GetCF(name, field, B_pars):
    '''Get CF energy levels'''
    ce = re(name, field, B_pars)
    return np.concatenate((np.array(ce.rawenergy[:, np.newaxis]), ce.moment), axis=1)


@xlfunc
def GetHParts(field):
    '''Get H direction and H size'''
    field = np.array(field)
    H_size = np.sqrt(np.dot(field, field.conj().transpose()))
    H_direction = field / H_size
    return np.reshape(H_direction.tolist()+[H_size,],(4,1))


@xlfunc
def GetSusc(energy, moment, H_direction, H_size, T):
    '''Get Susceptibility'''
    return susc(np.array(energy), np.array(moment), np.array(H_direction) , H_size, T)

@xlfunc
def CFConvert(fromUnit, toUnit):
    '''Returns conversion factor from fromUnit to toUnit.'''
    if fromUnit == "T/uB" and toUnit == "mol/emu":
        return C.C3
    if fromUnit == "mol/emu" and toUnit == "T/uB":
        return 1 / C.C3
