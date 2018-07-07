#!/usr/bin/env python
'''
test file
'''

import numpy as np
import pytest
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
    ising.build_ising(rand = False)
    b = ising.config
    assert(np.allclose(a, b))

def test_MC():
    L = 8
    temp = 2.0
    ising = Ising(L, temp) 
    ising.MC_kernel(init_steps = 5000, bin_steps = 100, mc_per_bin = 100, DEBUG = False)

if __name__ == '__main__':
    test_MC()



