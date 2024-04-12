import matplotlib.pyplot as plt
import numpy as np

from data import vitadata1
from logdelay import logdelay

parameters_1 = [5]  # adagolások száma / sűrűsége
parameters_2 = [0.5, 0.75]  # D_0 kezdeti mennyiség
for i in range(len(parameters_1)):
    for j in range(len(parameters_2)):
        t, r, y, g, m, d = logdelay([parameters_1[i], parameters_2[j]])

        fig, ax1 = plt.subplots(figsize=(8, 6))

        ax1.set_ylabel('Nr. of cells', color="black")
        ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
        ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
        ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
        ax1.tick_params(axis='y', labelcolor="black")
        ax1.set_xlabel('Time')

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-ax1
        ax2.set_ylabel('Drug level', color="blue")
        ax2.step(t, d, "blue", where='pre', alpha=0.4)
        ax2.tick_params(axis='y', labelcolor="blue")


        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
