import numpy as np
import matplotlib.pyplot as plt


def logdelay(params):
    treatment_freq, drug_0 = params
    r_p = 1/0.05
    k_r = 10
    k_y = 10
    k_g = 10
    cc_r = 1.028
    cc_y = 0.387
    cc_g = 1.285
    drug_decay = np.log(2) / 2  # derived from Doxorubicin half-life
    delta = 0.8  # drug mortality rate from Song 2022
    K_c_r = 0.25  # Let us try this as selectivity
    K_c_y = 0
    K_c_g = 0
    m0 = 0
    r0 = 106
    y0 = 209
    g0 = 76
    repeat = 20
    r_length, y_length, g_length = 10.28, 3.87, 12.85  # átlagos sejtciklus hossz becsülve a vitadello cikkben.
    K = 20000  # meg kellett emelni, mert elérte. (K - 2 * M * (len(R) + len(Y) + len(G))) < 0!!
    end_time = 47.75
    time_step = 0.25
    drug_step = 0.2
    Data = []
    #############################################

    for i in range(repeat):
        M = m0
        D = drug_0
        R = np.random.uniform(0, r_length, r0)
        Y = np.random.uniform(0, y_length, y0)
        G = np.random.uniform(0, g_length, g0)

        t = 0
        t_stat = [0]
        t_drug = [0]
        R_stat, Y_stat, G_stat, M_stat, D_stat = [len(R)], [len(Y)], [len(G)], [M], [drug_0]

        while t < end_time and np.any([0 < len(R), 0 < len(Y), 0 < len(G), (M + 2 * (len(R) + len(Y) + len(G))) < K]):

            # Gillespie
            a_r = r_p * M * (K - (M + 2 * (len(R) + len(Y) + len(G)))) / K
            a_d_r = K_c_r * len(R) * (1 - np.exp(-delta * D))
            a_d_y = K_c_y * len(Y) * (1 - np.exp(-delta * D))
            a_d_g = K_c_g * len(G) * (1 - np.exp(-delta * D))
            a = [a_r, a_d_r, a_d_y, a_d_g]
            a0 = sum(a)
            Tau = np.random.exponential(scale=(1 / a0))

            theta1 = R[np.argmin(R)] if 0 < len(R) else np.Inf
            theta2 = Y[np.argmin(Y)] if 0 < len(Y) else np.Inf
            theta3 = G[np.argmin(G)] if 0 < len(G) else np.Inf

            # Itt össze kell hasonlítani Tau-t a teták minimumával (azaz min(theta1, theta2, theta3)al)
            if Tau < min(theta1, theta2, theta3):
                min_tau = Tau

                # generate random number in U(0,1).
                r_1 = np.random.uniform(0, 1)
                # choose next reaction
                mu = 0  # this will be the index of the next reaction
                N = r_1 * a0 - a[mu]

                while N > 0:
                    mu = mu + 1
                    N = N - a[mu]

                if mu == 0:  # a_r nyert -> újabb sejt kerül a sejtciklusba
                    M -= 1
                    R -= min_tau
                    Y -= min_tau
                    G -= min_tau
                    time_in_R = np.random.gamma(k_r, cc_r)  # a sejtciklus R-ben töltött véletlen hossza
                    R = np.append(R, [time_in_R])
                elif mu == 1:  # a_d_r nyert
                    if len(R) > 0:
                        choice = np.random.choice(len(R), size=1)
                        R = np.delete(R, choice)
                elif mu == 2:  # a_d_y nyert
                    if len(Y) > 0:
                        choice = np.random.choice(len(Y), size=1)  # itt fogyott el Y
                        Y = np.delete(Y, choice)
                else:  # a_d_g nyert
                    if len(G) > 0:
                        choice = np.random.choice(len(G), size=1)
                        G = np.delete(G, choice)

            elif np.min(np.array([theta1, theta2, theta3])) == theta1:  # Ha a legrövidebb idő R-beli sejthez tartozik
                min_tau = theta1
                R = np.delete(R, np.argmin(R))
                time_in_Y = np.random.gamma(k_y, cc_y)  # a sejtciklus Y-ben töltött hossza (jöhet eloszlásból)
                Y = np.append(Y, [time_in_Y])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif np.min(np.array([theta1, theta2, theta3])) == theta2:  # Ha a legrövidebb idő Y-beli sejthez tartozik
                min_tau = theta2
                Y = np.delete(Y, np.argmin(Y))
                time_in_G = np.random.gamma(k_g, cc_g)  # a sejtciklus G-ben töltött hossza (jöhet eloszlásból)
                G = np.append(G, [time_in_G])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif np.min(np.array([theta1, theta2, theta3])) == theta3:  # Ha a legrövidebb idő G-beli sejthez tartozik
                min_tau = theta3
                G = np.delete(G, np.argmin(G))
                M += 2
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            t += min_tau
            drug_disc_step = int((t - t_drug[-1]) / drug_step)

            if drug_disc_step >= 0:
                for i in range(drug_disc_step):
                    t_drug.append(t_drug[-1] + drug_step)
                    if t_drug[-1] % treatment_freq == 0:
                        D += drug_0
                    else:
                        D = drug_0 * np.exp(-drug_decay * (t % treatment_freq))

            ######################## Statistics ################################

            disc_step = int((t - t_stat[-1]) / time_step)

            if disc_step >= 0:
                for j in range(disc_step):
                    t_stat.append(t_stat[-1] + time_step)
                    R_stat.append(len(R))
                    Y_stat.append(len(Y))
                    G_stat.append(len(G))
                    M_stat.append(M)
                    D_stat.append(D)

        ####################################################################
        Data.append([t_stat, R_stat, Y_stat, G_stat, M_stat, D_stat])

    m = len(Data[0][0])
    size = len(Data)
    for d in Data:
        if len(d[0]) < m:
            m = len(d[0])

    avg_t, avg_R, avg_Y, avg_G, avg_M, avg_D = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)

    for d in Data:
        avg_t += np.array(d[0])[:m] / size
        avg_R += np.array(d[1])[:m] / size
        avg_Y += np.array(d[2])[:m] / size
        avg_G += np.array(d[3])[:m] / size
        avg_M += np.array(d[4])[:m] / size
        avg_D += np.array(d[5])[:m] / size
    return [avg_t, avg_R, avg_Y, avg_G, avg_M, avg_D]

"""
# A szaporodási ráta inverze, a pihenő fázisban töltött idő várható értéke
Tau = .05
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
t, r, y, g, m, d = logdelay([1 / Tau, k_r, cc_r, k_y, cc_y, k_g, cc_g])

fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_ylabel('Nr. of cells', color="black")
ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
ax1.tick_params(axis='y', labelcolor="black")
ax1.set_xlabel('Time')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-ax1
ax2.set_ylabel('Drug level', color="blue")
ax2.step(t, d, "blue", where='post', alpha=0.4)
ax2.tick_params(axis='y', labelcolor="blue")


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
"""
