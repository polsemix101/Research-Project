import obspy
import ordpy
import numpy as np
import matplotlib.pyplot as plt

path = "C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Data/D/"

data = obspy.read(path+"2024.011.1717.53.GLOS.443.evt")

data.plot()

for j in range(0,3):
    data1 = data[j][:]
    n = len(data1)

    s = np.arange(0,n-200,200)
    e = np.arange(200,n,200)

    E = np.array([])
    for i in range(len(s)):
        E = np.append(E,ordpy.complexity_entropy(data1[s[i]:e[i]],dx=4,probs=False)[0])

    plt.plot(E)
    plt.show()