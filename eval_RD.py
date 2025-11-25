#!/usr/local/bin/python3

import sys
import subprocess
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

dist = list()
bitlen = list()
labels = list()
for predMode in range(0, 5):
    for quantSize in range(8, 65, 8):
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
plt.xlim(-50, None)
for simplex in simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
for vertex in vertices[0:8]:
    plt.text(points[vertex, 0] - 3, points[vertex, 1], str(labels[vertex]), horizontalalignment="right")
for qs in [16, 32, 64]:
    x_series = np.array(bitlen)[[labels[i]["quantSize"] == qs for i in range(len(points))]]
    y_series = np.array(dist)[[labels[i]["quantSize"] == qs for i in range(len(points))]]
    x_tickz = np.array([x_series.min() - 3, x_series.max() + 6])
    A = np.vstack([x_series, np.ones(len(x_series))]).T
    m, c = np.linalg.lstsq(A, y_series, rcond=None)[0]
    plt.plot(x_tickz, m * x_tickz + c, 'r:')
plt.colorbar()
plt.xlabel("Länge Bitstream [Kilobytes]")
plt.ylabel("Mittlerer quadratischer Fehler [1/Pixel]")

if (len(sys.argv) == 1):
    plt.show()
else:
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(sys.argv[1], format="svg", bbox_inches="tight", dpi=100)