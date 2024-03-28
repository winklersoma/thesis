# Az a kérdés, hogyan függjön a particio a bomlási rátától, illetve esetleg a drog szinttől.
# Érdemes lenne lépcsős függvényként ábrázolni először is és kipróbálni, hogy mik azok, amikre rábólintanánk.
# Másodsorban meg lehetne nézni a Peti által javasolt Linda JS Allen könyvben, hátha írnak valami okosságot.
# Ha mindez megvan és döntöttünk, akkor a kódban ki kell javítani az adagolást és a bomlás-updateket,
# Valamint a diszkretizálást is meg kell csinálni azzal a módszerel, ahogy a Peti által küldött Colabban szerpel.
import numpy as np
import matplotlib.pyplot as plt

drug_0 = 3  # Ez a mylogdelayből jön.
drug_decay = 0.9  # derived in 2022 Song (it is 0.5199 in dePillis 2009)
t_0 = 0
t_end = 5  # Feltettem, hogy pl. 5 óránként szeretnénk adagolni a gyógyit.
n_partitions = 21  # Próbaérték


partition = np.linspace(start=t_0, stop=t_end, num=n_partitions)
continuous = np.linspace(start=t_0, stop=t_end, num=1001)

out_continuous = drug_0 * np.exp(-drug_decay * continuous)

# Standard Deviation from the Mean
out_std_plus = out_continuous + np.sqrt(drug_0 * (np.exp(-drug_decay * continuous) - np.exp(-continuous)))
out_std_minus = out_continuous - np.sqrt(drug_0 * (np.exp(-drug_decay * continuous) - np.exp(-continuous)))

# Code for the staircase function
out_stairs = out_continuous.copy()
i = t_0 + 1
p_counter = t_end / (n_partitions-1)
while continuous[i] < t_end:
    if continuous[i] < p_counter:
        out_stairs[i] = out_stairs[i-1]
        i += 1
    else:
        out_stairs[i] = out_continuous[i]
        p_counter += (t_end / (n_partitions-1))
        i += 1


plt.plot(continuous, out_continuous, color='red', linestyle='dashed', label='Theoretical (expected) drug level')
plt.plot(continuous, out_stairs, color='green', label='Drug level in model')
plt.plot(continuous, out_std_plus, color='purple', label='STD of actual drug level')
plt.plot(continuous, out_std_minus, color='purple')

plt.title(' t_0=' + str(t_0) + ', t_end=' + str(t_end) +", number of partitions= " + str(n_partitions) + ", step size= " + str(t_end / (n_partitions - 1)))
plt.legend()
plt.xlabel("Time")
plt.ylabel("Drug level")
plt.show()
