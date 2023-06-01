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
plt.rcParams['figure.dpi']=180

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
            


P=float(Find('PB') or 0)                                  # days
Pdot=float(Find('PBDOT') or 0)                          #unitless
wdot=float(Find('OMDOT') or 0)                              # degrees/year
e=float(Find('ECC') or 0)                                #unitless
gamma= float(Find('GAMMA') or 0)                           #seconds 
x=float(Find('XDIST') or 0)                                  #light-seconds
s=float(Find('SINI') or 0)                                 #unitless
r=float(Find('SHAPMAX') or 0)                           #seconds*M_sun

print(float(Find('SINI') or 0),x )

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





fig, ax = plt.subplots()       
plt.xlim(0,3)
plt.ylim(0,3)
plt.xlabel("Pulsar mass")
plt.ylabel("Companion mass")

P1=ax.contour(mp, mc, r_(mp,mc)-r, [0],colors='orange')
P2=ax.contour(mp, mc, s_(mp,mc)-s, [0],colors='green')
P3=ax.contour(mp, mc, gamma_(mp,mc)-gamma, [0],colors='purple')    #contour plots of the curves
P4=ax.contour(mp, mc, P_dot(mp,mc)-Pdot, [0],colors='blue')
P5=ax.contour(mp, mc, w_dot(mp,mc)-wdot, [0],colors='red')

plt.title('Mass-mass diagram for ' + Find('PSRJ'))    #title with pulsar name 

h1,_ = P1.legend_elements()
h2,_ = P2.legend_elements()
h3,_ = P3.legend_elements()
h4,_ = P4.legend_elements()
h5,_ = P5.legend_elements()

ax.legend([h1[0], h2[0], h3[0], h4[0], h5[0]], ['Shapiro range', 'sin(i)','Einstein delay','Period dot',
           'omega dot']
)  
plt.show()