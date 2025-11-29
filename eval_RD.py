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
            #command = "./comp_demo radial1024.bin 1024 1024 " + str(predMode) + " " + str(quantSize) + " " + str(splitLevel)
            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            p.wait()
            output = output.decode("utf-8")
            dist_val = float(output.split("\n")[1].split(" ")[-1])
            bitlen_val = float(output.split("\n")[2].split(" ")[-1])
            #to make the figure looking nicer, remove higher bitrates
            if bitlen_val <= 5:
                dist.append(dist_val)
                bitlen.append(bitlen_val)
                labels.append({"predMode": predMode, "quantSize": quantSize, "splitLevel": splitLevel})

points = np.array([list(t) for t in zip(bitlen, dist)])
hull = ConvexHull(points)

simplices = list()
vertices = list()
for simplex in hull.simplices[1:]:
    p0 = [points[simplex, 0][0], points[simplex, 1][0]]
    p1 = [points[simplex, 0][1], points[simplex, 1][1]]
    slope = (p1[1] - p0[1]) / (p1[0] - p0[0])
    if (slope <= -10) and (slope > -100000):
        vertex0 = np.where(points == p0)[0][0]
        vertex1 = np.where(points == p1)[0][0]
        simplices.append(simplex)
        vertices = list(set(vertices + [vertex0, vertex1]))
assert(len(vertices) - 1 == len(simplices))
vertices = sorted(vertices, key=lambda x: labels[x]["quantSize"])

hull_parameters = {"qs": [], "lambda": [], "minCost": []}
for i, vertex in enumerate(vertices):
    qs = labels[vertex]["quantSize"]
    p = points[vertex]
    pl = None
    pr = None
    slope = None
    if (i > 0) and (i < len(vertices) - 1):
        pl = points[vertices[i - 1]]
        pr = points[vertices[i + 1]]
        slope = (pr[1] - pl[1]) / (pr[0] - pl[0])
    elif i == 0:
        pr = points[vertices[i + 1]]
        slope = (pr[1] - p[1]) / (pr[0] - p[0])
    elif i == len(vertices) - 1:
        pl = points[vertices[i - 1]]
        slope = (p[1] - pl[1]) / (p[0] - pl[0])
    minCost = dist[vertex] - slope * bitlen[vertex]
    hull_parameters["qs"].append(qs)
    hull_parameters["lambda"].append(-slope)
    hull_parameters["minCost"].append(minCost)
    print(str(points[vertex]) + ":  " + str(labels[vertex]))
    print("Lambda: " + str(-slope))
    print("Min. Costs: " + str(minCost))

costs = np.zeros(len(points))
for i in range(len(points)):
    qs = labels[i]["quantSize"]
    Lambda = hull_parameters["lambda"][hull_parameters["qs"].index(qs)]
    minJ = hull_parameters["minCost"][hull_parameters["qs"].index(qs)]
    costs[i] = (dist[i] + Lambda * bitlen[i]) / minJ

for figMode in range(2):
    cats = np.ones(len(points))
    title = ""
    if figMode == 0:
        cats = [t["quantSize"] for t in labels]
        title = "R-D-points for all compression modes grouped by quantization stepsize"
    elif figMode == 1:
        cats = costs
        title = "Increasing R-D-costs (J = D + λR) orthogonal to the optimal RD-curve"
    plt.figure()
    plt.scatter(bitlen, dist, c=cats, cmap="viridis_r")
    plt.xlim(-1.5, None)
    for simplex in simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    for vertex in vertices:
        plt.text(points[vertex, 0] - 0.1, points[vertex, 1], str(labels[vertex]), horizontalalignment="right")
    plt.colorbar()
    plt.xlabel("Mittlere Code-Länge [Bits/Pixel]")
    plt.ylabel("Mittlerer quadratischer Fehler [1/Pixel]")
    plt.title(title)

    if (len(sys.argv) == 1):
        plt.show()
    else:
        fig = plt.gcf()
        fig.set_size_inches(15, 9)
        fname = sys.argv[1]
        assert(fname.endswith(".svg"))
        if figMode == 1:
            fname = fname[:-4] + "_costs.svg"
        fig.savefig(fname, format="svg", bbox_inches="tight", dpi=100)