import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from protocol import logdelay

parameters_1 = [(0.9, 0, 0), (0, 0.9, 0), (0, 0, 0.9), (0.3, 0.3, 0.3), (0.45, 0, 0.45)]
parameters_2 = [1, 3, 5, 7, 9]

hm_data = np.zeros((len(parameters_1), len(parameters_2)))

for i in range(len(parameters_1)):
    for j in range(len(parameters_2)):
        t, r, y, g, m, d = logdelay(k_c=parameters_1[i], drug_0=parameters_2[j])
        total_unhealthy = r + y + g + m
        hm_data[i][j] = total_unhealthy[-1] / total_unhealthy[0]

df = pd.DataFrame(hm_data, index=parameters_1, columns=parameters_2)
fig, ax = plt.subplots(figsize=(8, 8))
hm = sns.heatmap(data=df, annot=True, cbar=False)
plt.xlabel('Dosage concentration in mg/l')
plt.ylabel('Sensitivity in compartments R, Y, G respectively')
plt.title("Relative count of unhealthy cells")
# plt.savefig('Sensitivity_vs_concentration', dpi=300)
plt.show()
