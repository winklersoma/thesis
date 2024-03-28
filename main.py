import matplotlib.pyplot as plt
import numpy as np

from data import vitadata1
from experimental import logdelay
# A szaporodási ráta inverze, a pihenő fázisban töltött idő várható értéke
Tau = .05
# A sejtciklus hossza most konstans, ha eloszlás lesz, azt máshol kell intézni és mindig újra sorsolni
Theta = 1
# A gamma eloszlások paraméterei becsülve a vitadello cikkből
# R kompartmentnél 10 alfázis esetén 10.28-as becsült fázishossz mellett k_r=10, cc_r=1.028:
k_r = 10
cc_r = 1.028
# Y kompartmentnél 10 alfázis esetén 3.87-es becsült fázishossz mellett k_y=10, cc_y=0.387:
k_y = 10
cc_y = 0.387
# G kompartmentnél 10 alfázis esetén 12.85-ös becsült fázishossz mellett k_g=10, cc_g=1.285:
k_g = 10
cc_g = 1.285
t, r, y, g, m = logdelay([1/Tau, k_r, cc_r, k_y, cc_y, k_g, cc_g])
plt.plot(t, r, c="red", linewidth=1.5)
plt.plot(t, y, c="orange", linewidth=1.5)
plt.plot(t, g, c="green", linewidth=1.5)
#plt.plot(t, m, c="brown", linewidth=1.5)
plt.plot(t, r+y+g+m, c="black", linewidth=1.5)
plt.scatter(t, vitadata1[0], c="red", alpha=0.15)
plt.scatter(np.arange(0, len(vitadata1[1])*0.25, 0.25), vitadata1[1], c="orange", alpha=0.15)
plt.scatter(np.arange(0, len(vitadata1[2])*0.25, 0.25), vitadata1[2], c="green", alpha=0.15)

plt.title("$r^{-1}= $" + str(Tau) + ", \t$\\vartheta= $" + str(Theta))
plt.ylabel("Nr. of cells")
plt.xlabel("Time")
# plt.savefig(str(Tau)+','+str(Theta)+'.png', dpi = 300)

plt.show()