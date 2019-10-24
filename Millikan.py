#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:43:02 2019
@author: espen
"""
import numpy as np
import uncertainties as u

#Defining Variables (all ufloat in (value, uncert.))
#--------------------
d = u.ufloat(5.47e-3, 0.02e-3) #Plate separation (m)
rho = u.ufloat(886, 0) #Oil density (kg*m^-3)
g = u.ufloat(9.8, 0) #Gravitational accel. (m*s^-2)
eta = u.ufloat(1.84e-5, 0) #Viscosity of air (N*s*m^-2)
b = u.ufloat(8.20e-3, 0) #Constant (Pa*m)
P = u.ufloat(1.01e5, 0.02e5)#Calculated pressure in lab (Pa)
h = u.ufloat(0.00050, 0.00001) #Distance between two marks (m)
V = u.ufloat(396.4, 2.0) #Voltage across plates (V)
#--------------------

#Importing rise/fall and finding mean +/- error
#--------------------
risefall = np.loadtxt("risefall.txt") #Load rising and fallind times (tab-delimited text)
rvals = risefall[:,0] #Col 1 is the rise time
fvals = risefall[:,1] #Col 2 is the fall time
rmean = np.mean(rvals) #Mean rise time
fmean = np.mean(fvals) #Mean fall time
rsd = np.std(rvals, ddof=1) #Rise standard deviation
fsd = np.std(fvals, ddof=1) #Fall standard deviation
riset = u.ufloat(rmean, rsd) #Mean rise time with error (s)
fallt = u.ufloat(fmean, fsd) #Mean fall time with error (s)
#--------------------
 
#Calculate rise, fall velocity and drop radius a
#--------------------
vr = h/riset #Rise velocity (m*s^-1)
vf = h/fallt #Fall velocity (m*s^-1)
print(vf)
print(vr)
a = ((9*eta*vf)/(2*g*rho))**(1/2) #Radius of drop (m)
#--------------------

#Calculate electron charge q in Coulombs 
#--------------------
#q = (4/3) * np.pi * d * 
q = (4/3)*np.pi*d*(((1/(g*rho))*(((9*eta)/2)**3))**(1/2))*((1/(1+(b/(P*a))))**(3/2))*(((vf + vr)*(vf**(1/2)))/V)
print('Radius a = {:P}'.format(a))
print('Charge q in C = {:P}'.format(q))

