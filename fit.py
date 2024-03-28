import numpy as np
from scipy.optimize import minimize

from mylogdelay import logdelay
from data import vitadata1

M0 = 0
R0 = vitadata1[0][0]
Y0 = vitadata1[1][0]
G0 = vitadata1[2][0]

# img2, WM983C cell line
vita_R = vitadata1[0]  # [:len(m)]  # erre azért van szükség, mert különben *-nál különböznek a numpy listák hosszai
vita_Y = vitadata1[1]  # [:len(m)]
vita_G = vitadata1[2]  # [:len(m)]  # és lehetne egy olyan konfiguráciot is kiprobálni, ahol máshogy vonok össze

# Miket nem akarok minimalizálni:
Repeat = 1  # szimuláció ismétlésszáma
Theta = 1  # A különböző kompartmentekben vett eredeti egyenletes eloszlás maximuma

# Miket akarok minimalizálni:
Tau = .05  # A szaporodási ráta inverze, a pihenő fázisban töltött idő várható értéke
r_p = 1 / Tau
'gamma eloszlások paraméterei'


def fn_min(params):
    t, r, y, g, m = logdelay(params)
    simu_R, simu_Y, simu_G = m + r, y, g
    return np.sum((simu_R - vita_R) ** 2) + np.sum((simu_Y - vita_Y) ** 2) + np.sum((simu_G - vita_G) ** 2)  # *


parameters = np.array([21.87278898, 31.03369373, 0.37221878, 31.94261499, 0.42267796, 30.84056807, 0.36664538])
print(fn_min(parameters))
res = minimize(fn_min, parameters, method='nelder-mead', options={'maxiter': 500})
print(res)
print(res.x)
