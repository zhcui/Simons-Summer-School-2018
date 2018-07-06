#!/usr/bin/env python
'''
Ising model for simons 

'''

import numpy as np
#import random as rd
from numpy import random as rd

### GLOBAL VARIABLE ###


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
        self.energy, self.mag = self.measure()

    def measure(self):
        A = self.config
        e = 0.0
        m = 0.0
	for x in xrange(self.length):
	    for y in xrange(self.length):
	        xl = x-1
		yl = y-1
		env_factor = A[x,y] * (A[xl,y] + A[x,yl])
		e -= env_factor
		m += A[x,y]
	return e, m
        
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
    

    def MC_kernel(self, init_steps = 5000, bin_steps = 100, mc_per_bin = 100):

        #e, m = self.measure()
        e = self.energy
        m = self.mag
        N = self.num_sites

        #init_steps=10000
        #bin_steps=100
        #mc_per_bin=100

        #pre equilibrium 

        for i in xrange(init_steps):
            self._MC_step()


        #measure

        E1 = 0.0
        E2 = 0.0
        M1 = 0.0
        M2 = 0.0

        for i in xrange(bin_steps):
            e1=0.0
            e2=0.0
            m1=0.0
            m2=0.0
            for j in xrange(mc_per_bin):
                self._MC_step()
                e1 += self.energy
                e2 += self.energy**2
                m1 += np.abs(self.mag)
                m2 += self.mag**2
            e1_ave = e1 / (float(mc_per_bin) * N)
            e2_ave = e2 / (float(mc_per_bin) * N**2)
            m1_ave = m1 / (float(mc_per_bin) * N)
            m2_ave = m2 / (float(mc_per_bin) *N**2)
            
            print [e1_ave,e2_ave,m1_ave,m2_ave]
            print self.config
            
            E1 += e1_ave
            E2 += e2_ave
            M1 += m1_ave
            M2 += m2_ave
                
        E1 /= float(bin_steps)
        E2 /= float(bin_steps)
        M1 /= float(bin_steps)
        M2 /= float(bin_steps)

        print "E1, E2, M1, M2"
        print E1, E2, M1, M2

        
    
