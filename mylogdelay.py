import numpy as np


def logdelay(params):
    r_p, k_r, cc_r, k_y, cc_y, k_g, cc_g = params
    drug_0 = 1000  # Why not
    drug_decay = 0.5199  # derived in dePillis 2009
    delta = 0.8  # drug mortality rate from Song 2022
    K_c = 0.1  # Let us try this as selectivity
    m0 = 0
    r0 = 106
    y0 = 209
    g0 = 76
    repeat = 1
    r_length, y_length, g_length = 10.28, 3.87, 12.85  # átlagos sejtciklus hossz becsülve a vitadello cikkben.
    K = 10000
    end_time = 47.75
    time_step = 0.25
    n_treatments = 10
    Data = []
    #############################################

    for i in range(repeat):
        M = m0
        D = drug_0
        R = np.random.uniform(0, r_length, r0)
        # Szerintem ez egy rossz döntés: ha több idejük lenne y-ban, és kevesebb G-ben, eredetileg, akkor jobban közelítene az igazi görbéhez.
        Y = np.random.uniform(0, y_length, y0)
        G = np.random.uniform(0, g_length, g0)

        discrete_time = time_step
        t = 0
        t_stat = [0]
        R_stat, Y_stat, G_stat, M_stat = [len(R)], [len(Y)], [len(G)], [M]

        while t < end_time and np.any([0 < len(R), 0 < len(Y), 0 < len(G), M < K]):

            tau1 = R[np.argmin(R)] if 0 < len(R) else np.Inf
            tau2 = Y[np.argmin(Y)] if 0 < len(Y) else np.Inf
            tau3 = G[np.argmin(G)] if 0 < len(G) else np.Inf
            tau4 = np.random.exponential(K / (M * r_p * (K - M - 2 * (len(R) + len(Y) + len(G))))) \
                if 0 < M * (K - M - 2 * (len(R) + len(Y) + len(G))) else np.Inf
            tau5 = K_c * (1 - np.random.exponential(scale=-delta * D, size=1)) * (len(R) + len(Y) + len(G))  # A sejtszám is önkényes.
            print(tau5, tau4)
            ####################
            min_tau = np.min(np.array([tau1, tau2, tau3, tau4]))
            if t < discrete_time:
                t += min_tau
            else:
                discrete_time += time_step
                D = drug_0 * np.exp(-drug_decay * t)
                t += min_tau
            print('Actual drug level = ' + str(D))

            if min_tau == tau1:  # Ha a legrövidebb idő R-beli sejthez tartozik
                R = np.delete(R, np.argmin(R))
                time_in_Y = np.random.gamma(k_y, cc_y)  # a sejtciklus Y-ben töltött hossza (jöhet eloszlásból)
                Y = np.append(Y, [time_in_Y])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau2:  # Ha a legrövidebb idő Y-beli sejthez tartozik
                Y = np.delete(Y, np.argmin(Y))
                time_in_G = np.random.gamma(k_g, cc_g)  # a sejtciklus G-ben töltött hossza (jöhet eloszlásból)
                G = np.append(G, [time_in_G])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau3:  # Ha a legrövidebb idő G-beli sejthez tartozik
                G = np.delete(G, np.argmin(G))
                M += 2
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau4:  # Enter cell cycle
                M -= 1
                R -= min_tau
                Y -= min_tau
                G -= min_tau
                time_in_R = np.random.gamma(k_r, cc_r)  # a sejtciklus hossza
                R = np.append(R, [time_in_R])

            ######################## Statistics ################################

            disc_step = int((t - t_stat[-1]) / time_step)

            if disc_step >= 0:
                for j in range(disc_step):
                    t_stat.append(t_stat[-1] + time_step)
                    R_stat.append(len(R))
                    Y_stat.append(len(Y))
                    G_stat.append(len(G))
                    M_stat.append(M)

        ####################################################################
        Data.append([t_stat, R_stat, Y_stat, G_stat, M_stat])

    m = len(Data[0][0])
    size = len(Data)
    for d in Data:
        if len(d[0]) < m:
            m = len(d[0])

    avg_t, avg_R, avg_Y, avg_G, avg_M = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)

    for d in Data:
        avg_t += np.array(d[0])[:m] / size
        avg_R += np.array(d[1])[:m] / size
        avg_Y += np.array(d[2])[:m] / size
        avg_G += np.array(d[3])[:m] / size
        avg_M += np.array(d[4])[:m] / size

    return [avg_t, avg_R, avg_Y, avg_G, avg_M]
