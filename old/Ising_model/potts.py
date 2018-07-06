#!/usr/bin/env python

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

### GLOBAL VARIABLE ###


L=10
N=L*L
temp=1.0

#pflip array
pflip=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
index=[-4,-2,0,2,4]
for i in index:
	pflip[i]=np.exp(-float(i)*2.0/(float(temp)))

pflip_potts=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
index_potts=[-4,-3,-2,-1,0,1,2,3,4]
for i in index_potts:
	pflip_potts[i]=np.exp(-float(i)*1.0/(float(temp)))

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
	

def build_potts(rand):
	if rand==1:
		Is=rd.randint(1,4,size=(L,L))
		return Is
	else:
		Is=rd.randint(1,2,size=(L,L))
		return Is
def measure_potts(A):
	e=0.0
	m=0.0
	M=0.0+0.0j
	for x in range(0,L):
		for y in range(0,L):
			xl=x-1
			yl=y-1
			curr=A[x,y]
			if A[xl,y]==curr:
				e=e-1.0
			if A[x,yl]==curr:
				e=e-1.0	
			m=m+float(curr)
			M=M+np.exp((curr-1)*2*np.pi*1j/3.0)
	return [e,m,M]

def MC_step_potts(A,e,m,M):
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
		curr=A[x,y]
		choice=rd.randint(1,3)
		if curr==1:
			trans=1+choice
		elif curr==2:
			if choice==1:
				trans=1
			else:
				trans=3
		else:
			trans=choice
		count_old=0
		count_new=0
		for j in [xl,xr]:
			if curr==A[j,y]:
				count_old=count_old+1
			elif trans==A[j,y]:
				count_new=count_new+1
		for j in [yl,yr]:
			if curr==A[x,j]:
				count_old=count_old+1
			elif trans==A[x,j]:
				count_new=count_new+1
		env_factor=count_old-count_new
		p=pflip_potts[env_factor]
		r=rd.random()
		if r<=p:
			A[x,y]=trans
			e=e+float(env_factor)
			m=m+float(trans-curr)
			M=M+np.exp((trans-1)*2*np.pi*1j/3.0)-np.exp((curr-1)*2*np.pi*1j/3.0)
	return [e,m,M]

### MAIN ###

	### variables ###
A=build_potts(1)
print A
[e,m,M]=measure_potts(A)
print [e,m,M]
print np.real(M)
print np.imag(M)

initsteps=1000
binstep=500
MC_per_bin=500

#pre eq

for i in xrange(0,initsteps):
	[e,m,M]=MC_step_potts(A,e,m,M)


#measure

E1=0.0
E2=0.0
M1=0.0
M2=0.0
MM=0.0+0.0j

for i in xrange(0,binstep):
	e1=0.0
	e2=0.0
	m1=0.0
	m2=0.0
	mm=0.0+0.0j
	for j in xrange(0,MC_per_bin):
		[e,m,M]=MC_step_potts(A,e,m,M)
		e1=e1+e
		e2=e2+e**2
		m1=m1+m
		m2=m2+m**2
		mm=mm+M
	e1_ave=e1/(float(MC_per_bin)*float(N))
	e2_ave=e2/(float(MC_per_bin)*float(N)**2)
	m1_ave=m1/(float(MC_per_bin)*float(N))
	m2_ave=m2/(float(MC_per_bin)*float(N)**2)
	mm_ave=mm/(float(MC_per_bin)*float(N))
	print [e1_ave,e2_ave,m1_ave,m2_ave]
	print mm_ave
	#print A
	E1=E1+e1_ave
	E2=E2+e2_ave
	M1=M1+m1_ave
	M2=M2+m2_ave
	MM=MM+mm_ave
	
E1=E1/float(binstep)
E2=E2/float(binstep)
M1=M1/float(binstep)
M2=M2/float(binstep)
MM=MM/float(binstep)

print [E1,E2,M1,M2]
print MM

Real=np.real(MM)
Imag=np.imag(MM)
print Real
print Imag
print "Order Parameter:"
print np.sqrt(Real**2+Imag**2)

'''
[e,m]=measure(A)



initsteps=10000
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
#	print A
	E1=E1+e1_ave
	E2=E2+e2_ave
	M1=M1+m1_ave
	M2=M2+m2_ave
	
E1=E1/float(binstep)
E2=E2/float(binstep)
M1=M1/float(binstep)
M2=M2/float(binstep)

print [E1,E2,M1,M2]

'''

#Is=build_ising(0)
#print Is
#Is=build_ising(0)
#print Is
#print pflip


