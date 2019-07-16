#!/usr/bin/env python3
from matplotlib import pyplot as plt
import numpy as np
x = np.linspace(0,10, num=100)
y = np.sin(x)
plt.plot(x,y, label="sin(x)")
plt.legend()
plt.show()