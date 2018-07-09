#!/usr/bin/env python

"""
2D Ising model with MC, for Simons Summer School 2018
Author: Zhihao Cui
        Gaurav Harsha 
"""

import numpy as np
#import random
from numpy import random as rd
import time

import os, sys
import importlib
import ctypes
from ctypes import *

libmc_tools_path = './libmc_tools.so'

if os.path.exists(libmc_tools_path):
    libmc_tools = CDLL('%s'%(libmc_tools_path)) 
else:
    print "libmc_tools.so not found."
    sys.exit(1)

c_mc_step = libmc_tools._Z9c_mc_stepPiRKiPKdRdS4_
c_mc_step_pre = libmc_tools._Z13c_mc_step_prePiRKiPKd


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

    def build_ising(self, rand = True, measure = True):
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

        if rand:
            # some are 1, some are -1
            A = rd.randint(-1, 1, size = (self.length, self.length), dtype = np.int32) # TODO: choose a more memory saving dtype
            A[A == 0] = 1
        else:
            # all one case
            A = rd.randint(1, 2, size = (self.length, self.length), dtype = np.int32)
        
        self.config = A
        
        if measure:
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
        
        e = 0.0
        m = 0.0
        for x in xrange(self.length):
            xl = x - 1
            for y in xrange(self.length):
                yl = y - 1
                env_factor = self.config[x, y] * (self.config[xl, y] + self.config[x, yl])
                e -= env_factor
                m += self.config[x, y]
        return e, m
    
    def _MC_step_pre(self):
        """
        One MC step for pre-equilibrium, only update the configuration.

        Parameters
        ----------

        Returns
        -------

        """
        A = self.config
            
        for i in xrange(self.num_sites):
            # randomly select a site
            x = rd.randint(0, self.length)
            y = rd.randint(0, self.length)
            
            # local interaction count, PBC is used
            env_factor = A[x, y] * (A[x - 1, y] + A[(x + 1)%self.length, y] + \
                         A[x, y - 1] + A[x, (y + 1)%self.length])
            # flip
            if rd.random() <= self.pflip[env_factor]:
                A[x, y] *= -1
    
    def _MC_step(self):
        """
        One MC step, update the configuration, energy and magnetization as well.

        Parameters
        ----------

        Returns
        -------

        """
        A = self.config
            
        for i in xrange(self.num_sites):
            # randomly select a site
            x = rd.randint(0, self.length)
            y = rd.randint(0, self.length)
            
            # local interaction count, PBC is used
            env_factor = A[x, y] * (A[x - 1, y] + A[(x + 1)%self.length, y] + \
                         A[x, y - 1] + A[x, (y + 1)%self.length])
            # flip
            if rd.random() <= self.pflip[env_factor]:
                A[x, y] *= -1
                self.energy += env_factor * 2.0
                self.mag += A[x, y] * 2.0

    def _MC_step_C(self):
        """
        One MC step, update the configuration, energy and magnetization as well.
        using C code;

        Parameters
        ----------

        Returns
        -------

        """
        e = c_double(self.energy)
        m = c_double(self.mag)
        c_mc_step(self.config.ctypes.data_as(c_void_p), byref(c_int(self.length)), self.pflip.ctypes.data_as(c_void_p), \
                  byref(e), byref(m))
        self.energy = e.value
        self.mag = m.value

    def _MC_step_pre_C(self):
        """
        One MC step, update the configuration, not energy and magnetization.
        using C code;

        Parameters
        ----------

        Returns
        -------

        """

        c_mc_step_pre(self.config.ctypes.data_as(c_void_p), byref(c_int(self.length)), self.pflip.ctypes.data_as(c_void_p))

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

        t0 = time.time()
        print "\n2D Ising model with MC algorithm\n"
        print "Parameters: length of lattice: %6d , temperature: %10.5f "%(self.length, self.temp)
        print "MC parameters: init_steps: %10d , bin_steps = %10d , mc_per_bin = %10d "\
               %(init_steps, bin_steps, mc_per_bin)

        if self.config is None:
            self.build_ising(rand = True, measure = False)

        N = self.num_sites

        #pre equilibrium 
        print "Pre-equilibrium..."
        t1 = time.time()

        for i in xrange(init_steps):
            self._MC_step_pre_C()
       
        self.energy, self.mag = self.measure()
    
        t2 = time.time()
        print "Pre-equilibrium finished, wall time %10.2f"%(t2 - t1)

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
                self._MC_step_C()
                e1_ave += self.energy
                e2_ave += self.energy**2
                m1_ave += np.abs(self.mag)
                m2_ave += self.mag**2
            e1_ave /= (float(mc_per_bin) * N)
            e2_ave /= (float(mc_per_bin) * N**2)
            m1_ave /= (float(mc_per_bin) * N)
            m2_ave /= (float(mc_per_bin) * N**2)
            
            print "  bin : %5d , <E> : %10.5f , <E^2> : %8.3f , <M> : %10.5f , <M^2> : %8.3f "\
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
        
        t3 = time.time()

        print "\nFinal results: \n"
        print "<E> : %12.6f , <E^2> : %10.5f , <M> : %12.6f , <M^2> : %10.5f "\
               %(E1, E2, M1, M2)
        print "total wall time: %15.3f s"%(t3 - t0)

        
    
