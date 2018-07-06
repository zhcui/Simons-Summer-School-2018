#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt

plt.figure(figsize=(9,6))
n=1000
x=np.random.randn(1,n)
y=np.random.randn(1,n)
T=np.arctan2(x,y)
plt.scatter(x,y,c=T,s=25,alpha=0.4,marker='o')
plt.show()
