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

from Ising import build_ising

def test_build_ising():
    L = 10
    a = np.ones((L, L))
    b = build_ising(L, rand = 0)
    assert(np.allclose(a, b))


