import numpy as np


def logdelay(params):
    r_p, cc_length, m0, p0, repeat = params

    K = 100000
    time_step = cc_length / 10
    Data = []
    #############################################

    for i in range(repeat):
        M = m0

        P = np.random.uniform(0, cc_length, p0)

        t = 0
        t_stat = [0]
        P_stat, M_stat = [len(P)], [M]

        while 0 < len(P) or M < K:  # len(P_stat) < 13: #

            if 0 < len(P):  # osztódó sejtek
                cell = np.argmin(P, axis=None, out=None)  # visszaadja a legkisebb elem indexét
                tau1 = P[cell]
            else:
                tau1 = np.Inf

            if 0 < M * (K - M - 2 * len(P)):
                tau2 = np.random.exponential(K / (M * r_p * (K - M - 2 * len(P))))
            else:
                tau2 = np.Inf

            ####################

            if tau1 < tau2:  # Division
                P = np.delete(P, cell)
                M += 2
                if len(P) != 0:
                    P -= tau1
                t += tau1

            elif tau2 <= tau1:  # Enter cell cycle
                M -= 1
                if len(P) != 0:
                    P -= tau2
                cell_cycle = cc_length  # a sejtciklus hossza, lehet pl. np.random.gamma(shape, scale)
                P = np.append(P, [cell_cycle])
                t += tau2

            ######################## Statistics ################################

            disc_step = int((t - t_stat[-1]) / time_step)

            if disc_step >= 0:
                for j in range(disc_step):
                    t_stat.append(t_stat[-1] + time_step)
                    P_stat.append(len(P))
                    M_stat.append(M)

        ####################################################################

        Data.append([t_stat, P_stat, M_stat])

    m = len(Data[0][0])
    size = len(Data)
    for d in Data:
        if len(d[0]) < m:
            m = len(d[0])

    avg_t, avg_P, avg_M = np.zeros(m), np.zeros(m), np.zeros(m)

    for d in Data:
        avg_t += np.array(d[0])[:m] / size
        avg_P += np.array(d[1])[:m] / size
        avg_M += np.array(d[2])[:m] / size

    return [avg_t, avg_P, avg_M]
