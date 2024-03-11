import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = np.loadtxt("rect/data/outflow_temp_1.csv")

temp = data[:]



dof_data = pd.read_csv('rect/data/outflow_temp_dofs.csv', index_col=0)
t = np.float64(dof_data.columns)
dofs = np.int64(dof_data.index)
temps_dof= dof_data.values

plt.plot(t, temp, label = 'reference solution')

for i, d in enumerate(dofs):
    plt.plot(t, temps_dof[i], label=f"{d} spatial dofs")

plt.legend()
plt.show()

for i, d in enumerate(dofs):
    plt.plot(t, np.abs(temps_dof[i]-temp), label=f"{d} spatial dofs")

plt.legend()
plt.yscale('log')
plt.show()

plt.plot(dofs, np.abs(temps_dof[:, -1]-temp[-1]), '.',  label=f"{d} spatial dofs")

plt.legend()
plt.yscale('log')
plt.show()