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
plt.rcParams['figure.dpi']=200  #sets the resoloution 

Ms=1.989e30
c=3e8
G=6.6743e-11


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def Find(param,n):       # fucntion that searches listed file for value of parameter
                                              # n gives the nth number after the parameter name
    with open(sys.argv[1], 'r') as fp:
        lines = fp.readlines()
        for line in lines:    
            if line.find(param) != -1:
                l= lines[lines.index(line)]
                l=l.strip() 
                l=l.split(" ")
                k=0
                for i in np.arange(0,len(l),1):                    
                    if l[i-k] =="":
                        del l[i-k]
                        k=k+1
                for i in np.arange(0,len(l),1):
                     if i>=n:
                         if is_number(l[i]):
                           return float(l[i])
                         else: return l[i]
        print('Could not find  ' + param)
        return 0 
    
# pulling parameters from file , first element is the value and second is the error

P=Find('PB',1)                                                        # days
Pdot=[Find('PBDOT',1),Find('PBDOT',2)]                          #unitless
wdot=[Find('OMDOT',1) ,Find('OMDOT',2)]                      # degrees/year
e=Find('ECC',1)                                           #unitless
gamma= [Find('GAMMA',1),Find('GAMMA',2)]                      #seconds 
x=Find('A1',1)                                               #light-seconds
s=[Find('SINI',1) ,Find('SINI',2)]                               #unitless
r=[Find('SHAPMAX',1) ,Find('SHAPMAX',2) ]                          #seconds*M_sun

print(P,Pdot,wdot,e,gamma,x,s,r)


if e==0:
    print('Warning: eccentricity cannot be 0')

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


'''
Mp=[]
Mc=[]
S=[]
R=[]
for i in np.arange(0.01,3,0.01):
 Mp.append(i)
 Mc.append(1)
 S.append(w_dot(i,1))       #stuff for debugging
 R.append(P_dot(i,1))
plt.plot(Mp,S)
plt.figure()
plt.plot(Mp,R)'''

fig, ax = plt.subplots()       
plt.xlim(0,3)
plt.ylim(0,3)
plt.xlabel("Pulsar mass")
plt.ylabel("Companion mass") 

#contour plots of the curves, width of the curves is accounted for by the error of the measured parameters
# If the error is 0 just a curve is plotted



if r[1]!= 0:
 P1=ax.contourf(mp, mc, r_(mp,mc)-r[0], [-r[1],r[1]],colors='orange',alpha=0.5)
else:
    P1=ax.contour(mp, mc, r_(mp,mc)-r[0], [0],colors='orange',alpha=0.5)
if s[1] != 0:
 P2=ax.contourf(mp, mc, s_(mp,mc)-s[0], [-s[1],s[1]],colors='green',alpha=0.5)
else:
    P2=ax.contour(mp, mc, s_(mp,mc)-s[0], [0],colors='green',alpha=0.5)
if gamma[1] != 0:
 P3=ax.contourf(mp, mc, gamma_(mp,mc)-gamma[0], [-gamma[1],gamma[1]],colors='purple',alpha=0.5)   
else:
   P3=ax.contour(mp, mc, gamma_(mp,mc)-gamma[0], [0],colors='purple',alpha=0.5)
if Pdot[1] != 0: 
 P4=ax.contourf(mp, mc, P_dot(mp,mc)-Pdot[0], [-Pdot[1],Pdot[1]],colors='blue',alpha=0.5)
else:
   P4=ax.contour(mp, mc, P_dot(mp,mc)-Pdot[0], [0],colors='blue',alpha=0.5)
if wdot[1] != 0:
 P5=ax.contourf(mp, mc, w_dot(mp,mc)-wdot[0], [-wdot[1],wdot[1]],colors='red',alpha=0.5)
else:
    P5=ax.contour(mp, mc, w_dot(mp,mc)-wdot[0], [0],colors='red',alpha=0.5)


plt.title('Mass-mass diagram for ' + Find('PSRJ',1))    #title with pulsar name 

h1,_ = P1.legend_elements()
h2,_ = P2.legend_elements()
h3,_ = P3.legend_elements()
h4,_ = P4.legend_elements()
h5,_ = P5.legend_elements()

ax.legend([h1[0], h2[0], h3[0], h4[0], h5[0]], ['Shapiro range', 'sin(i)','Einstein delay','Period dot',
           'omega dot'])

plt.show()