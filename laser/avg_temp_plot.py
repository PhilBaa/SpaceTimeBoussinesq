import numpy as np
import matplotlib.pyplot as plt

fine = np.genfromtxt('laser/data/avg_temp.csv')

plt.plot(fine)
plt.show()

for i in range(5):
    d = np.genfromtxt(f'laser/data/avg_temp_{i}.csv')
    plt.plot(np.abs(fine - d))
plt.show()