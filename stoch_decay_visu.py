import numpy as np
import matplotlib.pyplot as plt
from stoch_exp_decay import exp_stoch

r = 0.9
D0 = 1000
dt = .02
time_step = 0.001
Data = []
for i in range(50):  # The simulation will run 100 times!
    d = exp_stoch(D0, r, time_step)
    Data.append(d)  # Elements of Data are list of lists
    # (Data is a list of lists of lists)

# minimal length
m = np.inf
for i in range(len(Data)):
    if len(Data[i][2]) < m: m = len(Data[i][2])  # Let m be the shortest time_stat length
print("m = " + str(m))

# make stat
stat_time = np.arange(0, time_step * (m), time_step)  # This will be the common time axis
print("stat_time = " + str(stat_time))
print(len(stat_time))

Data_mtx = []  # Each row of Data_mtx will consist of the states corresponding to the discrete times
for i in range(len(Data)):
    Data_mtx.append(Data[i][3][:m])
print(np.shape(Data_mtx))
print(5 * "_")
Data_mtx = np.array(Data_mtx)
mean = []
var = []
for i in range(np.shape(Data_mtx)[1]):
    mean.append(np.mean(Data_mtx[:, i]))
    var.append(np.var(Data_mtx[:, i]))

mean = np.array(mean)
var = np.array(var)

print(len(mean))

print(5 * "--")

DT1 = np.arange(0, time_step * (m), dt)  # DT1 is a new time axes with dt-long intervals

print(DT1)
DT = []
DE = []
for i in range(len(DT1) - 1):  # Staircase function code!
    DT.append(DT1[i])
    DT.append(DT1[i + 1])
    DE.append(D0 * np.exp(-r * DT1[i]))
    DE.append(D0 * np.exp(-r * DT1[i]))

for i in range(len(Data)):
    plt.plot(Data[i][0], Data[i][1], alpha=0.5)

plt.title('Exponential decay of ' + str(1000) + ' drug molecules at a decay rate of ' + str(r))
plt.plot(stat_time, mean, linewidth=2, color="purple", label='Mean drug level')
# plt.plot(stat_time, D0*np.exp(-4*stat_time), linewidth = 1, color = "white")
# plt.fill_between(stat_time, mean - np.sqrt(var), mean + np.sqrt(var), color='b', alpha=.95)
plt.plot(stat_time, mean + np.sqrt(var), linewidth = 2.5, color="black", label='Mean level plus standard deviation')
plt.plot(stat_time, mean - np.sqrt(var), linewidth = 1.5, color="black")
plt.plot(DT, DE, linewidth=1.2, color="red", label='Expected drug level in model')
plt.xlabel("Time")
plt.ylabel("Drug level")
plt.legend()
plt.show()