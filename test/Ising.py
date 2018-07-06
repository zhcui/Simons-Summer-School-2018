#!/usr/bin/env python
'''
Ising model for simons 

'''


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


#L = 10
#N = L * L
#temp = 1.5
#
##pflip array
#pflip=np.zeros(9, dtype = np.double)
#index=[-4,-2,0,2,4]
#for i in index:
#    pflip[i] = np.exp(-i*2.0/temp)
#

class Ising(object):
    def __init__(self, length, temp): # TODO multi dim
        self.length = length
        self.num_sites = length * length
        self.temp = temp
        #pflip array
        self.pflip = np.zeros(9, dtype = np.double)
        index=[-4, -2, 0, 2, 4]
        if temp <= 1e-10:
            raise ValueError('Temperature should be greater than zero, should be larger than 1e-10',temp)
        for i in index:
            self.pflip[i] = np.exp(-i * 2.0 / temp)

        self.energy = 0.0
        self.mag = 0.0
        self.config = None

    def build_ising(self, rand):
        """
        Build the Ising model in a 2D array.

        Parameters
        ----------
        rand : bool
            whether randomly generate or not

        Returns
        -------
        np.ndarray
            The Ising model configurations.

        """
        L = self.length

        if rand == 1:
            # some are 1, some are -1
            Is = rd.randint(-1, 1, size = (L, L)) # TODO: choose a more memory saving dtype
            Is[Is==0] = 1
        else:
            # all one case
            Is=rd.randint(1,2,size=(L,L))
        self.config = Is
        
    def _MC_step(self):
            
        e = self.energy
        m = self.mag
        A = self.config
        L = self.length

        for i in xrange(0, self.num_sites):
            # randomly select a site
            x=rd.randint(0,L)
            y=rd.randint(0,L)
           
            # boundary PBC
            xl = np.mod(x-1,L)
            xr = np.mod(x+1,L)
            yl = np.mod(y-1,L)
            yr = np.mod(y+1,L)
            # local interaction count
            env_factor = A[x,y]*(A[xl,y]+A[xr,y]+A[x,yl]+A[x,yr])
            p = self.pflip[env_factor]
            r = rd.random()
            # flip
            if r<=p:
                A[x,y] *= -1
                e += env_factor * 2.0
                m += A[x,y] * 2.0

        self.energy = e
        self.mag = m
        self.config = A
    






### FUNCTION ###

def build_ising(L, rand):
    """
    Build the Ising model in a 2D array.

    Parameters
    ----------
    rand : bool
        whether randomly generate or not

    Returns
    -------
    np.ndarray
        The Ising model configurations.

    """

    if rand == 1:
        # some are 1, some are -1
	Is = rd.randint(-1,1,size=(L,L)) # TODO: choose a more memory saving dtype
        Is[Is==0] = 1
	return Is
    else:
        # all one case
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
		env_factor = A[x,y] * (A[xl,y] + A[xr,y] + A[x,yl] + A[x,yr])
		p = self.pflip[env_factor]
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

if __name__ == '__main__':


    ### variables ###
    A=build_ising(1)
    print A

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



