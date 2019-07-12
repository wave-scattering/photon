from matplotlib import pylab as plt
import numpy as np
filename = "rtaAgcyl.dat"

data = np.loadtxt(filename)
print(data)
plt.scatter(data[:,0],data[:,1])
plt.show()
