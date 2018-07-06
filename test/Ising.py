#!/usr/bin/env python

"""
2D Ising model with MC, for Simons Summer School 2018
Author: Zhihao Cui
        Gaurav Harsha 
"""

import numpy as np
#import random as rd
from numpy import random as rd


class Ising(object):
    def __init__(self, length, temp): # TODO multi dim
        self.length = length
        self.num_sites = length * length
        self.temp = temp
        #pflip array
        self.pflip = np.zeros(9, dtype = np.double)
        index=[-4, -2, 0, 2, 4]
        if temp <= 1e-10:
            raise ValueError('Temperature should be greater than zero, should be larger than 1e-10', temp)
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
            Whether randomly generate or not

        Returns
        -------
        np.ndarray
            The Ising model configurations.

        """
        L = self.length

        if rand:
            # some are 1, some are -1
            Is = rd.randint(-1, 1, size = (L, L)) # TODO: choose a more memory saving dtype
            Is[Is == 0] = 1
        else:
            # all one case
            Is = rd.randint(1, 2, size = (L, L))
        self.config = Is
        self.energy, self.mag = self.measure()

    def measure(self):
        """
        Calculate the energy and magnetization of a given Ising configuration.

        Parameters
        ----------

        Returns
        -------
        double, double
            The energy and magnetization.

        """
        A = self.config
        e = 0.0
        m = 0.0
	for x in xrange(self.length):
	    for y in xrange(self.length):
	        xl = x - 1
		yl = y - 1
		env_factor = A[x, y] * (A[xl, y] + A[x, yl])
		e -= env_factor
		m += A[x, y]
	return e, m
        
    def _MC_step(self):
        """
        One MC step, update the configuration, energy and magnetization as well.

        Parameters
        ----------

        Returns
        -------

        """
            
        e = self.energy
        m = self.mag
        A = self.config
        L = self.length

        for i in xrange(self.num_sites):
            # randomly select a site
            x = rd.randint(0, L)
            y = rd.randint(0, L)
            # boundary PBC
            xl = np.mod(x - 1, L)
            xr = np.mod(x + 1, L)
            yl = np.mod(y - 1, L)
            yr = np.mod(y + 1, L)
            # local interaction count
            env_factor = A[x, y] * (A[xl, y] + A[xr, y] + A[x, yl] + A[x, yr])
            p = self.pflip[env_factor]
            r = rd.random()
            # flip
            if r<=p:
                A[x, y] *= -1
                e += env_factor * 2.0
                m += A[x, y] * 2.0

        self.energy = e
        self.mag = m
        self.config = A
    
    def MC_kernel(self, init_steps = 5000, bin_steps = 100, mc_per_bin = 100, DEBUG = False):
        """
        Main MC procedure, measure <E>, <E^2>, <M>, <M^2> during the MC.

        Parameters
        ----------
        init_steps : int
            Number of steps in pre-equilibrium procedure.
        bin_steps : int
            Number of bins for measurement.
        mc_per_bin : int
            Number of MC steps per bin.
        DEBUG: bool
            Whether to print out DEBUG information.

        Returns
        -------

        """

        print "\n2D Ising model with MC Algorithm\n"
        print "Parameters: length of lattice: %6d , temperature: %10.5f "%(self.length, self.temp)
        print "MC parameters: init_steps: %10d , bin_steps = %10d , mc_per_bin = %10d "\
            %(init_steps, bin_steps, mc_per_bin)

        e = self.energy
        m = self.mag
        N = self.num_sites

        #pre equilibrium 
        print "Pre-equilibrium..."
        for i in xrange(init_steps):
            self._MC_step()
        print "Pre-equilibrium finished"

        #measure
        E1 = 0.0
        E2 = 0.0
        M1 = 0.0
        M2 = 0.0
        
        # loop of bins
        for i in xrange(bin_steps):
            e1_ave = 0.0
            e2_ave = 0.0
            m1_ave = 0.0
            m2_ave = 0.0

            # loop of mc per bin
            for j in xrange(mc_per_bin):
                self._MC_step()
                e1_ave += self.energy
                e2_ave += self.energy**2
                m1_ave += np.abs(self.mag)
                m2_ave += self.mag**2
            e1_ave /= (float(mc_per_bin) * N)
            e2_ave /= (float(mc_per_bin) * N**2)
            m1_ave /= (float(mc_per_bin) * N)
            m2_ave /= (float(mc_per_bin) * N**2)
            
            print "bin : %5d , <E> : %10.5f , <E^2> : %8.3f , <M> : %10.5f , <M^2> : %8.3f "\
                  %(i, e1_ave, e2_ave, m1_ave, m2_ave)
            if DEBUG:
                print self.config
                #e, m = self.measure()
                #assert((e == self.energy) and (m == self.mag))
            
            E1 += e1_ave
            E2 += e2_ave
            M1 += m1_ave
            M2 += m2_ave
                
        E1 /= float(bin_steps)
        E2 /= float(bin_steps)
        M1 /= float(bin_steps)
        M2 /= float(bin_steps)

        print "\nFinal Results: \n"
        print "<E> : %12.6f , <E^2> : %10.5f , <M> : %12.6f , <M^2> : %10.5f "\
              %(E1, E2, M1, M2)

        
    
