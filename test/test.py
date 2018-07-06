#!/usr/bin/env python
'''
test file
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

import pytest

#from Ising import build_ising
from Ising import Ising

def test_init_ising_zero_temp():
    L = 10
    temp = 0.0
    with pytest.raises(ValueError):
        ising = Ising(L, temp) 

def test_build_ising():
    L = 10
    temp = 1.0
    ising = Ising(L, temp) 
    a = np.ones((L, L))
    b = ising.build_ising(rand = 0)
    assert(np.allclose(a, b))



