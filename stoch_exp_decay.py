import numpy as np

"""
Simulates an exponential decay from D0 with rate r and
returns with a list of lists consisting of the
[sampled cumulative random times,
the states in the random times,
the discrete times,
and the corresponding states (no. of molecules) ].
"""


def exp_stoch(D0, r, time_step):
    D = D0

    state = [D]
    time = [0]

    statistic = [D]
    time_stat = [0]
    while time_stat[-1] < 10 + time_step:  # D > 0:
        a = r * D
        Tau = np.random.exponential(scale=(1 / a))
        D -= 1
        state.append(D)
        time.append(time[-1] + Tau)

        disc_step = int((time[-1] - time_stat[-1]) / time_step)

        if disc_step > 0:
            for j in range(disc_step):
                time_stat.append(time_stat[-1] + time_step)
                statistic.append(D)

    return [time, state, time_stat, statistic]