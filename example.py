from crysfipy.reion import re, susceptibility as susc
import crysfipy.const as C
import numpy as np


print("*** CrysFiPy test ***")
print()
print("Cerium 3+")
ce = re("Ho", [1,0,0], 
        ["t", -0.173477508,
        0.001084591,
        -0.012701252,
        -3.34835E-06,
        0.0000097,
        ])
ce.calculate()
ce.getlevels()
print("{0.J2p1:.0f} energy levels, J = {0.J}, gJ = {0.gJ}\nalpha = {0.Alpha}, beta = {0.Beta}, gamma = {0.Gamma}".format(C.ion("ce")))
print()
print(ce)
print()
print("Calculation of susceptibility")
temps = [5,10,50,100,300]
for T in temps:
    print("T = {0} K \tchi_CF = {1} uB/T".format(T, susc(ce, T)))