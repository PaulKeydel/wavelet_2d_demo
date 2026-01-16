#!/usr/local/bin/python3

import sys
import subprocess
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

def collect_data() -> tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    dist = list()
    bitlen = list()
    qs = list()
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
                    qs.append(quantSize)
                    labels.append({"predMode": predMode, "quantSize": quantSize, "splitLevel": splitLevel})
    return np.array(bitlen), np.array(dist), np.array(qs), labels

def get_conv_hull(points: np.ndarray, qs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
    vertices = sorted(vertices, key=lambda x: qs[x])
    return np.array(vertices), np.array(simplices)

def calc_slopes(points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    slopes = np.empty(len(points))
    slopes.fill(np.nan)
    for i, vertex in enumerate(vertices):
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
        slopes[vertex] = slope
        minCost = points[vertex, 1] - slope * points[vertex, 0]
        print("V" + str(i) + " : " + str(points[vertex]))
        print("  Lambda: " + str(-slope))
        print("  Min. Costs: " + str(minCost))
    return slopes

def interpolate_lambda(vertices: np.ndarray, qs: np.ndarray, slopes: np.ndarray) -> np.ndarray:
    qs_hull = list(qs[vertices])
    lambda_hull = list(-(slopes[vertices]))
    z = np.polyfit(qs_hull, lambda_hull, 2)
    quad_fit = np.poly1d(z)
    lambda_vec = np.vectorize(lambda t: quad_fit(t))
    print("Lambda prediction from quantization stepsize:")
    print(quad_fit)
    print("Predicted Lambda values:")
    print([quad_fit(t) for t in qs_hull])
    return lambda_vec(qs)

def plot_RD(save_as: str | None):
    bitlen, dist, qs, labels = collect_data()
    points = np.column_stack((bitlen, dist))
    vertices, simplices = get_conv_hull(points, qs)
    slopes = calc_slopes(points, vertices)
    lambdas = interpolate_lambda(vertices, qs, slopes)
    costs = dist + lambdas * bitlen
    minJ = np.array([costs[qs == t].min() for t in qs])
    costs /= minJ

    for figMode in range(2):
        cats = np.ones(len(points))
        title = ""
        if figMode == 0:
            cats = qs
            title = "R-D-points for all compression modes grouped by quantization stepsize"
        elif figMode == 1:
            cats = costs
            title = "Increasing R-D-costs (J = D + λR) orthogonal to the optimal RD-curve"
        fig = plt.figure(figsize=(10, 6))
        plt.scatter(bitlen, dist, c=cats, cmap="viridis_r")
        plt.xlim(-1.8, None)
        for simplex in simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], "k-")
        for vertex in vertices:
            plt.text(points[vertex, 0] - 0.1, points[vertex, 1], str(labels[vertex]), horizontalalignment="right")
            if figMode == 0:
                continue
            p = points[vertex]
            m_orth = -1 / slopes[vertex]
            x_tickz = np.array([p[0] - 0.1, p[0] + 1.0])
            plt.plot(x_tickz, m_orth * (x_tickz - p[0]) + p[1], "g:")
        #plt.gca().set_yscale("log")
        plt.colorbar()
        plt.xlabel("Mittlere Code-Länge [Bits/Pixel]")
        plt.ylabel("Mittlerer quadratischer Fehler [1/Pixel]")
        plt.title(title)

        if (save_as is None):
            plt.tight_layout()
            plt.show()
        else:
            fig.set_size_inches(15, 9)
            assert(save_as.endswith(".svg"))
            fname = save_as if figMode == 0 else save_as[:-4] + "_costs.svg"
            fig.savefig(fname, format="svg", bbox_inches="tight", dpi=100)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        plot_RD(None)
    if len(sys.argv) == 2:
        plot_RD(sys.argv[1])