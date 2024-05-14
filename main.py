import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from logdelay import logdelay

# Define the parameters to run the model with
parameters_1 = [24, 48]  # Treatment frequency
parameters_2 = [15, 30]  # D_0 dosage size

# Create variable names for subplots
var_names = []
for i in range(len(parameters_1)):
    for j in range(len(parameters_2)):
        var_names.append('sub' + str(i) + str(j))
print(var_names)

# Initialize output
hm_data = np.zeros((len(parameters_1), len(parameters_2)))
fig = plt.figure(figsize=(9, 9))

k = 1
for i in range(len(parameters_1)):  # treatment frequency in rows
    for j in range(len(parameters_2)):  # dosage size in cols
        t, r, y, g, m, d = logdelay(treatment_freq=parameters_1[i], drug_0=parameters_2[j])
        total_unhealty = r + y + g + m
        hm_data[i][j] = total_unhealty[-1] / total_unhealty[0]

        var_names[k-1] = plt.subplot(len(parameters_1), len(parameters_2), k)
        var_names[k-1].plot(t, total_unhealty)
        k += 1
plt.suptitle("Trajectories of total nr. of cells with different parametrization")
plt.show()
df = pd.DataFrame(hm_data, index=parameters_1, columns=parameters_2)
hm = sns.heatmap(data=df, annot=True)
plt.xlabel('Dosage size (in ng)')
plt.ylabel('Treatment frequency (hours between two treatments)')
plt.title("Relative count of unhealty cells")
plt.show()


fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_ylabel('Nr. of cells', color="black")
ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
ax1.plot(t, r+y+g, c='black', linewidth=2, label='Total nr. of cells')
ax1.tick_params(axis='y', labelcolor="black")
ax1.set_xlabel('Time')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-ax1
ax2.set_ylabel('Drug level', color="blue")
ax2.step(t, d, "blue", where='pre', alpha=0.4)
ax2.tick_params(axis='y', labelcolor="blue")
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
