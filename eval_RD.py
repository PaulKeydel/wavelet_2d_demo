#!/usr/local/bin/python3

import os
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
        self.range_quant = np.arange(4, 25, 4)
        self.range_lagr  = np.arange(25, 1550, 75)
        for quantSize in self.range_quant:
            for lagrMult in self.range_lagr:
                os.chdir("comp_demo")
                command = "./comp_demo ../" + binImg + " " + str(width) + " " + str(height) + " " + str(quantSize) + " " + str(lagrMult)
                p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                (output, err) = p.communicate()
                p.wait()
                os.chdir("..")
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
        idx_min_rate = np.argmin(self.bitlen)
        idx_min_dist = np.argmin(self.dist)
        pr = self.points[idx_min_rate]
        pd = self.points[idx_min_dist]
        lin_thres = np.vectorize(lambda rate: (pr[1] - pd[1]) * (rate - pr[0]) / (pr[0] - pd[0]) + pr[1])

        fmask = self.dist < lin_thres(self.bitlen)
        fmask[idx_min_rate] = True
        fmask[idx_min_dist] = True
        idx_map = np.arange(len(self.points))[fmask]

        hull = ConvexHull(self.points[fmask])
        vertices = idx_map[hull.vertices]
        vertices = np.array(sorted(vertices, key=lambda x: self.bitlen[x]))
        simplices = np.array([idx_map[s] for s in hull.simplices])

        idx_ex = np.where((simplices == [idx_min_rate, idx_min_dist]).all(axis=-1) | (simplices == [idx_min_dist, idx_min_rate]).all(axis=-1))[0]
        simplices = np.delete(simplices, idx_ex, axis=0)
        assert(len(vertices) - 1 == len(simplices))
        return vertices, simplices

    def calc_slopes(self, vertices: np.ndarray) -> np.ndarray:
        range_qs  = self.range_quant
        slopes    = np.empty(len(self.points))
        qs_hull   = self.qs[vertices]
        dist_hull = self.dist[vertices]
        rate_hull = self.bitlen[vertices]

        slopes.fill(np.nan)
        for i, vertex in enumerate(vertices):
            qs = qs_hull[i]
            qs_prev = qs
            qs_next = qs
            if qs > range_qs.min():
                qs_prev = range_qs[np.where(range_qs == qs)[0][0] - 1]
            if qs < range_qs.max():
                qs_next = range_qs[np.where(range_qs == qs)[0][0] + 1]
            rate0 = rate_hull[qs_hull == qs_prev].mean()
            dist0 = dist_hull[qs_hull == qs_prev].mean()
            rate1 = rate_hull[qs_hull == qs_next].mean()
            dist1 = dist_hull[qs_hull == qs_next].mean()
            slopes[vertex] = (dist1 - dist0) / (rate1 - rate0)

        return slopes

    def interpolate_lambda(self, vertices: np.ndarray, slopes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        qs_hull     = self.qs[vertices]
        lambda_hull = -slopes[vertices]

        z = np.polyfit(qs_hull, lambda_hull, 2)
        quad_fit = np.poly1d(z)
        lambdas = np.apply_along_axis(lambda t: quad_fit(t), 0, self.qs)

        costs = self.dist + lambdas * self.bitlen
        minJ = np.array([costs[self.qs == t].min() for t in self.qs])
        costs /= minJ
        return lambdas, costs


if __name__ == "__main__":
    rd = RDeval("astronaut.bin", 512, 512)
    vertices, simplices = rd.get_conv_hull()
    slopes = rd.calc_slopes(vertices)
    lambdas, costs = rd.interpolate_lambda(vertices, slopes)

    d = {"rate": rd.bitlen, "dist": rd.dist, "qs": rd.qs, "slopes": slopes, "lambdas": -slopes, "lambda_pred": lambdas, "costs": costs}
    df_full = pd.DataFrame(data=d)
    df_hull = df_full.iloc[vertices]

    print(df_full.drop(["slopes", "lambdas"], axis=1).sort_values(["qs", "dist"], ascending=[True, True]).to_string())
    print()

    print(df_hull.to_string())
    print()
    print("Lambda prediction from quantization stepsize:")
    print(np.poly1d(np.polyfit(df_hull["qs"], df_hull["lambdas"], 2)))