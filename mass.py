# -*- coding: utf-8 -*-
"""
Created on Mon May 29 12:04:57 2023

@author: thoma
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import arange
from numpy import meshgrid
import sys
plt.rcParams['figure.dpi']=200

Ms=1.989e30
c=3e8
G=6.6743e-11


def Find(param):       # fucntion that searches file for value of parameter

    with open(sys.argv[1], 'r') as fp:
        lines = fp.readlines()
        for line in lines:    
            if line.find(param) != -1:
                l= lines[lines.index(line)]
                l=l.strip()
                l=l.split(" ")
                for i in np.arange(0,len(l),1):
                    if l[i] !="" and i>=1:
                        return l[i]
            
print(Find('F0'))

P=0.19765096149    # days
Pdot=-3.8e-13      #unitless
wdot=5.3105        # degrees/year
e=0.171875         #unitless
gamma= 7.9e-4      #seconds 
x=1.85891          #light-seconds
s=(2*.82)/(.82**2+1) #unitless
r=(2.7e-6)/(.82**3)  #seconds*M_sun

Ts=G*Ms/c**3
f=(1+(73/24)*e**2+(37/96)*e**4)/((1-e**2)**(7/2))


def P_dot(mp,mc):
    
    k=((-192*np.pi)/5)*((Ts*2*np.pi/(P*(60*60*24)))**(5/3))*f
    return k*(mp*mc/((mp+mc)**(1/3)))
   
def w_dot(mp,mc):
    
    k = 3*(Ts)**(2/3)*(2*np.pi/(P*(60*60*24)))**(5/3)/(1-e**2)
    k=k*365*60*60*24*180/np.pi
    return k*(mp+mc)**(2/3)

def gamma_(mp,mc):
    
    k=((Ts)**(2/3))*(((P*(60*60*24))/(2*np.pi))**(1/3))*e
    
    return k*mc*(mp+2*mc)/((mp+mc)**(4/3))

def s_(mp,mc):
    
    k=(Ts)**(-1/3)*(((P*(60*60*24))/(2*np.pi))**(-2/3))*x
    
    return k*((mp+mc)**(2/3))/mc

def r_(mp,mc):
     return Ts*mc




delta = 0.001
xrange = arange(0.01, 3.0, delta)
yrange = arange(0.01, 3.0, delta)
mp, mc = meshgrid(xrange,yrange)


Mp=[]
W=[]
Mc=[]
for i in np.arange(0.01,3,0.01):
    Mp.append(i)
    Mc.append(1)
    W.append(w_dot(i,1))


#plt.plot(Mp,W)
plt.figure()        

plt.contour(mp, mc, r_(mp,mc)-r, [0],colors='orange')
plt.contour(mp, mc, s_(mp,mc)-s, [0],colors='green')
plt.contour(mp, mc, gamma_(mp,mc)-gamma, [0],colors='purple')
plt.contour(mp, mc, P_dot(mp,mc)-Pdot, [0],colors='blue')
plt.contour(mp, mc, w_dot(mp,mc)-wdot, [0],colors='red')
plt.show()