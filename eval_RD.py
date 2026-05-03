#!/usr/local/bin/python3

import subprocess
from scipy.spatial import ConvexHull
import numpy as np
import pandas as pd

class RDeval:
    def __init__(self, binImg: str, width: int, height: int):
        dist    = list()
        bitlen  = list()
        qs      = list()
        lambdas = list()
        for quantSize in range(4, 33, 4):
            for lagrMult in range(50, 2501, 50):
                command = "./comp_demo " + binImg + " " + str(width) + " " + str(height) + " " + str(quantSize) + " " + str(lagrMult)
                p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                (output, err) = p.communicate()
                p.wait()
                output = output.decode("utf-8")
                dist.append(float(output.split("\n")[1].split(" ")[-1]))
                bitlen.append(float(output.split("\n")[2].split(" ")[-1]))
                qs.append(quantSize)
                lambdas.append(lagrMult)
        self.bitlen = np.array(bitlen)
        self.dist   = np.array(dist)
        self.qs     = np.array(qs)
        self.lambdas = np.array(lambdas)
        self.points = np.column_stack((bitlen, dist))

    def get_conv_hull(self) -> tuple[np.ndarray, np.ndarray]:
        hull = ConvexHull(self.points)
        vertices = hull.vertices.copy()
        vertices = np.roll(vertices, -np.argmin(self.bitlen[vertices]))
        vertices = vertices[:(np.argmin(self.dist[vertices]) + 1)] #not exclude R-D point with lowest distortion from convex hull

        simplices = list()
        for simplex in hull.simplices:
            p0 = [self.points[simplex, 0][0], self.points[simplex, 1][0]]
            p1 = [self.points[simplex, 0][1], self.points[simplex, 1][1]]
            if (p0 in self.points[vertices]) and (p1 in self.points[vertices]):
                simplices.append(simplex)
        assert(len(vertices) - 1 == len(simplices))
        #vertices = sorted(vertices, key=lambda x: self.qs[x])
        return vertices, np.array(simplices)

    def calc_slopes(self, vertices: np.ndarray) -> np.ndarray:
        slopes = np.empty(len(self.points))
        slopes.fill(np.nan)
        for i, vertex in enumerate(vertices):
            p = self.points[vertex]
            pl = None
            pr = None
            slope = None
            if (i > 0) and (i < len(vertices) - 1):
                pl = self.points[vertices[i - 1]]
                pr = self.points[vertices[i + 1]]
                slope = (pr[1] - pl[1]) / (pr[0] - pl[0])
            elif i == 0:
                pr = self.points[vertices[i + 1]]
                slope = (pr[1] - p[1]) / (pr[0] - p[0])
            elif i == len(vertices) - 1:
                pl = self.points[vertices[i - 1]]
                slope = (p[1] - pl[1]) / (p[0] - pl[0])
            slopes[vertex] = slope
        return slopes

    def interpolate_lambda(self, vertices: np.ndarray, slopes: np.ndarray) -> tuple[np.ndarray, np.poly1d]:
        qs_hull = list(self.qs[vertices])
        lambda_hull = list(-(slopes[vertices]))
        z = np.polyfit(qs_hull, lambda_hull, 2)
        quad_fit = np.poly1d(z)
        lambda_vec = np.vectorize(lambda t: quad_fit(t))
        return lambda_vec(self.qs), quad_fit

    def calc_costs(self, lambdas: np.ndarray) -> np.ndarray:
        costs = self.dist + lambdas * self.bitlen
        minJ = np.array([costs[self.qs == t].min() for t in self.qs])
        costs /= minJ
        return costs


if __name__ == "__main__":
    rd = RDeval("astronaut.bin", 512, 512)
    vertices, simplices = rd.get_conv_hull()
    slopes = rd.calc_slopes(vertices)
    lambdas, quad_fit = rd.interpolate_lambda(vertices, slopes)
    costs = rd.calc_costs(lambdas)

    d = {"rate": rd.bitlen, "dist": rd.dist, "qs": rd.qs, "slopes": slopes, "lambdas": -slopes, "lambdas_qs": lambdas, "costs": costs}
    df = pd.DataFrame(data=d)

    print(df.drop(["slopes", "lambdas"], axis=1).sort_values(["qs", "dist"], ascending=[True, True]).to_string())
    print()
    print(df.iloc[vertices].to_string())
    print()
    print("Lambda prediction from quantization stepsize:")
    print(quad_fit)