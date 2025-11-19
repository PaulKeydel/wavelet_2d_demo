#!/usr/local/bin/python3

import subprocess
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

dist = list()
bitlen = list()
labels = list()
for predMode in range(0, 5):
    for quantSize in range(4, 65, 4):
        for splitLevel in range(2, 6):
            command = "./comp_demo astronaut.bin 512 512 " + str(predMode) + " " + str(quantSize) + " " + str(splitLevel)
            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            p.wait()
            output = output.decode("utf-8")
            #print(output)
            dist.append(float(output.split("\n")[1].split(" ")[-1]))
            bitlen.append(int(output.split("\n")[2].split(" ")[-1]) // 1024)
            labels.append({"predMode": predMode, "quantSize": quantSize, "splitLevel": splitLevel})

points = np.array([list(t) for t in zip(bitlen, dist)])
hull = ConvexHull(points)
simplices = hull.simplices[3:]
vertices = hull.vertices[3:]

print("Points forming the Pareto frontier:")
for vertex in vertices:
    print(str(points[vertex]) + ":  " + str(labels[vertex]))

cats = [t["quantSize"] for t in labels]
plt.scatter(bitlen, dist, c=cats)
for simplex in simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.colorbar()
plt.xlabel("Länge Bitstream [Kilobytes]")
plt.ylabel("Mittlerer quadratischer Fehler [1/Pixel]")
plt.show()