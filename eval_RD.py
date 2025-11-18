#!/usr/local/bin/python3

import subprocess
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

dist = list()
bitlen = list()
qsize = list()
for predMode in range(0, 5):
    for quantSize in range(4, 65, 4):
        for splitLevel in range(2, 6):
            command = "./comp_demo astronaut.bin 512 512 " + str(predMode) + " " + str(quantSize) + " " + str(splitLevel)
            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            p.wait()
            output = output.decode("utf-8")
            #print(output)
            qsize.append(quantSize)
            dist.append(float(output.split("\n")[1].split(" ")[-1]))
            bitlen.append(int(output.split("\n")[2].split(" ")[-1]) // 1024)
#print(dist)
#print(bitlen)

points = np.array([list(t) for t in zip(bitlen, dist)])
hull = ConvexHull(points)

plt.scatter(bitlen, dist, c=qsize)
for simplex in hull.simplices[3:]:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.colorbar()
plt.xlabel("Länge Bitstream [Kilobytes]")
plt.ylabel("Mittlerer quadratischer Fehler [1/Pixel]")
plt.show()