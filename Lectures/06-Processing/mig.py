"""Zero-offset migration"""
import sys
import re
import struct
import numpy as np
import matplotlib.pyplot as pl
from math import *

#-----------------------------
# functions
#-----------------------------
        
def ricker(fp,tp,M,dt) :
   s=np.zeros((M)) 
   for k in range(0,M):
       t=k*dt
       wp=2.0*3.14159*fp
       s[k] =(1-0.5*pow(wp,2)*pow((t-tp),2))*exp(-0.25*pow(wp,2)*pow(t-tp,2)) 
   return(s)

#---------------
# Main script
#--------------

#Create a model

Nx=81
Nz=101
Nr=Nx
dt=0.004
Nt=500
Ns=50
dx=12.5
dz=12.5
L=(Nx-1)*dx
D=(Nz-1)*dz
z0=500.0
d=400.0
PI=3.14159
c0=2000.0
model=np.zeros((Nx))
data=np.zeros((Nr,Nt))
pulse=ricker(30.0,0.1,Nt,dt)
xa=np.zeros(Nx)
#===========================
#Modeling
#===========================

for i in range(0,Nx) :
    x=dx*i
    xa[i]=x
    a=(PI)/L
    model[i] = z0+d*sin(a*x) 

for i in range(0,Nx) :
    for j in range(0,Nr): 
        z=model[i]
        l=sqrt(pow(z,2)+pow((i-j)*dx,2))
        t = 2*l/c0
        it = int(t/dt)
        for k in range(0,Ns):
            if(it+k < Nt) :
                data[j,it+k]=data[j,it+k]+pulse[k]


print "*****Start mig"
#===========================
#Migration
#===========================
za=np.zeros(Nz)
img=np.zeros((Nx,Nz))
for i in range(0,Nx) :
    x=dx*i
    print "x : ", x
    for j in range(0,Nz): 
        z=dx*j
        za[i]=x
        for k in range(0,Nr) :
            l=sqrt(pow(z,2)+pow((i-k)*dx,2))
            t = 2*l/c0
            it = int(t/dt)
            if(it < Nt) :
                img[i,j]=data[k,it]+img[i,j]

#================
#Plot figures
#================
print np.amax(data)
fig=pl.figure()
pl.plot(xa,model)
pl.xlim(0,1000)
pl.ylim(0,1100)
pl.gca().invert_yaxis()
pl.xlabel("Distance (m)")
pl.ylabel("Depth (m)")
pl.savefig('model.pdf')
pl.show()
#pl.plot(pulse)
im=pl.imshow(np.transpose(np.flip(data,axis=1)),clim=(-2.5,2.5),extent=(0.0,1000.0,1.0,0.0),cmap='gray')
#Set aspect ratio
ar=1.0
ax=pl.gca()
asr = 1.0/(ax.get_data_ratio()*ar)
pl.Axes.set_aspect(ax,asr)
pl.xlabel("Distance (m)")
pl.ylabel("Time (sec)")
pl.savefig('data.pdf')
pl.show()

im=pl.imshow(np.transpose(img),extent=(0,1000.0,1000.0,0.0),cmap='gray')
ar=1.0
ax=pl.gca()
asr = 1.0/(ax.get_data_ratio()*ar)
pl.Axes.set_aspect(ax,asr)
pl.xlabel("Distance (m)")
pl.ylabel("Depth (m)")
pl.savefig('mig.pdf')
pl.show()

