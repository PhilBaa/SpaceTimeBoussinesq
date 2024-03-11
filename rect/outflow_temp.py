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

plt.xlabel('time $t$')
plt.ylabel(r'$J(\theta(t))$')
plt.xlim((0, t[-1]))
plt.ylim((0, 0.09))

plt.legend()
plt.savefig('rect/plots/outflow_temp_dofs.png')
plt.close()

fig, ax = plt.subplots(1, 2, figsize=(15, 8))

for i, d in enumerate(dofs):
    ax[0].plot(t, np.abs(temps_dof[i]-temp)/temp, label=f"{d} spatial dofs")

ax[0].legend()
ax[0].set_yscale('log')
ax[0].set_xlim((0, t[-1]))
ax[0].set_xlabel('time $t$')
ax[0].set_ylabel('rel. temperature error')
ax[0].set_title('Error of QoI with respect to time')

ax[1].plot(dofs, np.abs(temps_dof[:, -1]-temp[-1]), '.')

ax[1].set_yscale('log')
ax[1].set_xlabel('spatial dofs')
ax[1].set_ylabel('abs. temperature error')
ax[1].set_title('Error of QoI at final time')
fig.savefig('rect/plots/outflow_temp_dofs_error.png')