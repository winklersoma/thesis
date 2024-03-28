import numpy as np
import matplotlib.pyplot as plt


scale, size = 2, 10000  # mean=4, std=2*sqrt(2)
s = np.random.exponential(scale=scale, size=size)
s = np.sort(s)
t = np.linspace(0, 1000, 10000)
print(s.mean())
plt.plot(t, s, linewidth=1, color='r')
plt.show()

# print(np.random.gamma(40.53336864, 9.25942445 / 40.53336864))
