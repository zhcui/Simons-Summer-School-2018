#!/usr/bin/env python
#Fsqt_v3.py    A script for F_s(q,t) calculation, should be combined with get_pos.sh and LAMMPS
#              You should specify the # of particles, lattice constant, first peak of RDF, begin and end index for calculation and the temp you want calculate.
#              get_pos.py is to treat the LAMMPS dummp file to get the posistion informations. You should specify the lines you want to hold.
#              This script will generate the F_s(q,t)-vs-logt data in plot_fsqt file.

import os
import numpy as np
import scipy
import string
import math
#import random as rd
import cmath
from math import floor
from numpy import linspace
from numpy import random as rd
import time

### GLOBAL VARIABLE ###


t0 = time.time()
L=8
N=L*L
temp = 2.0

#pflip array
pflip=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
index=[-4,-2,0,2,4]
for i in index:
	pflip[i]=np.exp(-float(i)*2.0/(float(temp)))


### FUNCTION ###

def build_ising(rand):
	if rand==1:
		Is=rd.randint(-1,1,size=(L,L))
		for i in range(0,L):
			for j in range(0,L):
				if Is[i,j]==0:
					Is[i,j]=1
		return Is
	else:
		Is=rd.randint(1,2,size=(L,L))
		return Is

def MC_step(A,e,m):
	for i in xrange(0,N):
		x=rd.randint(0,L)
		y=rd.randint(0,L)
		if x==L-1:
			xl=L-2
			xr=0
		else:
			xl=x-1
			xr=x+1
		if y==L-1:
			yl=L-2
			yr=0
		else:
			yl=y-1
			yr=y+1
		env_factor=A[x,y]*(A[xl,y]+A[xr,y]+A[x,yl]+A[x,yr])
		p=pflip[env_factor]
		r=rd.random()
		if r<=p:
			A[x,y]=-A[x,y]
			e=e+float(env_factor*2)
			m=m+float(A[x,y]*2)
	return [e,m]

def measure(A):
	e=0.0
	m=0.0
	for x in range(0,L):
		for y in range(0,L):
			xl=x-1
			yl=y-1
			env_factor=A[x,y]*(A[xl,y]+A[x,yl])
			e=e-float(env_factor)
			m=m+float(A[x,y])
	return [e,m]
	



### MAIN ###

	### variables ###
def main():
    A=build_ising(1)
    print A

    [e,m]=measure(A)



    initsteps=5000
    binstep=100
    MC_per_bin=100

    #pre eq

    for i in xrange(0,initsteps):
            [e,m]=MC_step(A,e,m)


    #measure

    E1=0.0
    E2=0.0
    M1=0.0
    M2=0.0

    for i in xrange(0,binstep):
            e1=0.0
            e2=0.0
            m1=0.0
            m2=0.0
            for j in xrange(0,MC_per_bin):
                    [e,m]=MC_step(A,e,m)
                    e1=e1+e
                    e2=e2+e**2
                    m1=m1+abs(m)
                    m2=m2+m**2
            e1_ave=e1/(float(MC_per_bin)*float(N))
            e2_ave=e2/(float(MC_per_bin)*float(N)**2)
            m1_ave=m1/(float(MC_per_bin)*float(N))
            m2_ave=m2/(float(MC_per_bin)*float(N)**2)
            print [e1_ave,e2_ave,m1_ave,m2_ave]
            print A
            E1=E1+e1_ave
            E2=E2+e2_ave
            M1=M1+m1_ave
            M2=M2+m2_ave
            
    E1=E1/float(binstep)
    E2=E2/float(binstep)
    M1=M1/float(binstep)
    M2=M2/float(binstep)

    print [E1,E2,M1,M2]



    t1 = time.time()
    print t1-t0


if __name__ == '__main__':
    main()
