import os
import subprocess
from pathlib import Path
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import pywt
from pylab import cm

def run_compression(binImg: str, width: int, height: int, predMode: int, quantSize: int, splitLevel: int) -> tuple[float, float]:
    command = "./comp_demo " + binImg + " " + str(width) + " " + str(height) + " " + str(predMode) + " " + str(quantSize) + " " + str(splitLevel)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    p.wait()
    output = output.decode("utf-8")
    out_dist = float(output.split("\n")[1].split(" ")[-1])
    out_bitlen = float(output.split("\n")[2].split(" ")[-1])
    return out_dist, out_bitlen

def load_binaries(width: int, height: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    num_pixel = width * height
    #read all binary files
    buffer = Path("orig.bin").read_bytes()
    orig = np.reshape(list(struct.unpack('i'*num_pixel, buffer)), (width, height))

    buffer = Path("reco.bin").read_bytes()
    reco = np.reshape(list(struct.unpack('i'*num_pixel, buffer)), (width, height))

    buffer = Path("pred.bin").read_bytes()
    pred = np.reshape(list(struct.unpack('i'*num_pixel, buffer)), (width, height))

    buffer = Path("coeffs.bin").read_bytes()
    coeff = np.reshape(list(struct.unpack('i'*num_pixel, buffer)), (width, height))

    buffer = Path("resi.bin").read_bytes()
    resi = np.reshape(list(struct.unpack('i'*num_pixel, buffer)), (width, height))

    return orig, reco, pred, resi, coeff

def load_encoding(fname: str, width: int, height: int) -> np.ndarray:
    lines = []
    with open(fname, "r", encoding="UTF-8") as txtfile:
        lines = txtfile.read().splitlines()
    assert(len(lines) == width * height)
    return np.reshape(np.array(lines, dtype=object), (width, height))

def calc_entropy(message: np.ndarray):
    msg = message.flatten(order='C')
    n_labels = len(msg)
    if n_labels <= 1:
        return 0
    value, counts = np.unique(msg, return_counts=True)
    probs = counts / n_labels
    n_classes = len(value)
    if n_classes <= 1:
        return 0
    ent = 0.
    for i in probs:
        ent -= i * math.log(i, 2)
    return ent


class DemoPrediction:
    @staticmethod
    def visualize(save_as: str = ""):
        #close current opened plot window
        plt.close()

        width = 7
        height = 7
        ref = np.array([(21, 96, t) for t in range(16, 256, 16)], dtype=int)
        img = np.ones((height + 2, width + 2, 3), dtype=int) * 255
        img[0,2:,:] = ref[-width:]
        img[2:,0,:] = ref[:height]

        fig, axs = plt.subplots(1, 3, figsize=(12, 8))

        img[2:,2:,:] = np.repeat(ref[:height][:,None,:], width, axis=1)
        axs[0].imshow(img)
        axs[0].set_title("Prediction 1: Horizontal")
        axs[0].xaxis.set_ticks([])
        axs[0].yaxis.set_ticks([])

        img[2:,2:,:] = np.repeat(ref[-width:][None,:,:], height, axis=0)
        axs[1].imshow(img)
        axs[1].set_title("Prediction 2: Vertical")
        axs[1].xaxis.set_ticks([])
        axs[1].yaxis.set_ticks([])

        for y in range(height):
            for x in range(width):
                img[2+y,2+x,:] = (ref[x-width] + ref[height-y-1]) // 2
        axs[2].imshow(img)
        axs[2].set_title("Prediction 3: Diagonal (mean)")
        axs[2].xaxis.set_ticks([])
        axs[2].yaxis.set_ticks([])

        if save_as == "":
            plt.show()
        else:
            assert(save_as.endswith(".svg"))
            print("Processing file '" + save_as + "'...")
            fig.savefig(save_as, format="svg", bbox_inches="tight")


class DemoSteps:
    @staticmethod
    def visualize(binImg: str, width: int, height: int, predMode: int, quantSize: int, splitLevel: int, save_as: str = ""):
        assert(splitLevel > 0)
        #close current opened plot window
        plt.close()

        run_compression(binImg, width, height, predMode, quantSize, splitLevel)
        orig, reco, pred, resi, coeff = load_binaries(width, height)

        print("  orig: max=" + str(np.max(orig)) + " min=" + str(np.min(orig)))
        print("  reco: max=" + str(np.max(reco)) + " min=" + str(np.min(reco)))
        print("  coef: max=" + str(np.max(coeff)) + " min=" + str(np.min(coeff)))
        print("  resi: max=" + str(np.max(resi)) + " min=" + str(np.min(resi)))

        entr_orig = calc_entropy(orig)
        entr_coeff = calc_entropy(coeff)

        gridsize = 2 ** (8 - splitLevel)
        clip_x0 = 0
        clip_x1 = 256
        clip_y0 = 0
        clip_y1 = 256
        fig, axs = plt.subplots(2, 3, figsize=(12, 8))

        img1 = orig
        axs[0, 0].matshow(img1, cmap=cm.gray)
        axs[0, 0].set_title("Original")
        axs[0, 0].xaxis.set_ticks_position("bottom")

        img2 = reco
        axs[1, 0].matshow(img2, cmap=cm.gray)
        axs[1, 0].set_title("(5) Reconstructed")
        axs[1, 0].xaxis.set_ticks([])
        axs[1, 0].yaxis.set_ticks([])

        img3 = orig[clip_y0:clip_y1, clip_x0:clip_x1]
        axs[0, 1].matshow(img3, cmap=cm.gray)
        axs[0, 1].set_title("(1) Partitioning")
        axs[0, 1].xaxis.set_ticks([])
        axs[0, 1].yaxis.set_ticks([])
        if clip_x0 != 0 or clip_x1 != width or clip_y0 != 0 or clip_y1 != height:
            for y in range(clip_y0, 1 + (clip_y1 - clip_y0) // 2, gridsize):
                axs[0, 1].plot([clip_x0, (clip_x1 - clip_x0) // 2], [y, y], linewidth=1, color='r')
            for x in range(clip_x0, 1 + (clip_x1 - clip_x0) // 2, gridsize):
                axs[0, 1].plot([x, x], [clip_y0, (clip_y1 - clip_y0) // 2], linewidth=1, color='r')

        img4 = np.abs(coeff[clip_y0:((clip_y0 + clip_y1) // 2), clip_x0:((clip_x0 + clip_x1) // 2)])
        axs[1, 1].matshow(img4, cmap=cm.gray_r)
        axs[1, 1].set_title("(4) Transformed and quantized")
        axs[1, 1].xaxis.set_ticks([])
        axs[1, 1].yaxis.set_ticks([])

        img5 = pred[clip_y0:((clip_y0 + clip_y1) // 2), clip_x0:((clip_x0 + clip_x1) // 2)]
        axs[0, 2].matshow(img5, cmap=cm.gray)
        axs[0, 2].set_title("(2) Prediction signal")
        axs[0, 2].xaxis.set_ticks([])
        axs[0, 2].yaxis.set_ticks([])

        img6 = np.abs(resi[clip_y0:((clip_y0 + clip_y1) // 2), clip_x0:((clip_x0 + clip_x1) // 2)])
        axs[1, 2].matshow(img6, cmap=cm.gray_r)
        axs[1, 2].set_title("(3) Residual (absolute values)")
        axs[1, 2].xaxis.set_ticks([])
        axs[1, 2].yaxis.set_ticks([])

        txt = "Entropy original image: " + "{:.3f}".format(entr_orig) + "   /   entropy transformed image: " + "{:.3f}".format(entr_coeff)
        plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

        if save_as == "":
            plt.show()
        else:
            assert(save_as.endswith(".svg"))
            print("Processing file '" + save_as + "'...")
            fig.savefig(save_as, format="svg", bbox_inches="tight")

class DemoEncoding:
    @staticmethod
    def visualize(binImg: str, width: int, height: int, predMode: int, quantSize: int, splitLevel: int, save_as: str = ""):
        #close current opened plot window
        plt.close()

        run_compression(binImg, width, height, predMode, quantSize, splitLevel)
        orig, reco, pred, resi, coeff = load_binaries(width, height)

        enc_orig = load_encoding("enc_orig.txt", width, height)
        enc_coeff = load_encoding("enc_comp.txt", width, height)

        #shape of the zoom window
        wsx = 60
        wsy = 60
        verge = 5
        totalx = 3 * verge + 2 * wsx
        totaly = 3 * verge + 2 * wsy
        maxtrafo = np.max(np.abs(coeff))
        fig, axs = plt.subplots(1, 2, figsize=(12, 8))

        x0 = 0
        y0 = 0
        img1 = 255 * (orig - orig.min()) / np.ptp(orig)
        img1[y0:(y0 + totaly), x0:(x0 + totalx)] = 255 * np.ones((totaly, totalx), dtype=int)
        img1[(y0 + verge):(y0 + verge + wsy), (x0 + verge):(x0 + verge + wsx)] = orig[y0, x0] * np.ones((wsy, wsx), dtype=int)
        img1[(y0 + 2 * verge + wsy):(y0 + 2 * verge + 2 * wsy), (x0 + verge):(x0 + verge + wsx)] = orig[y0 + 1, x0] * np.ones((wsy, wsx), dtype=int)
        img1[(y0 + verge):(y0 + verge + wsy), (x0 + 2 * verge + wsx):(x0 + 2 * verge + 2 * wsx)] = orig[y0, x0 + 1] * np.ones((wsy, wsx), dtype=int)
        img1[(y0 + 2 * verge + wsy):(y0 + 2 * verge + 2 * wsy), (x0 + 2 * verge + wsx):(x0 + 2 * verge + 2 * wsx)] = orig[y0 + 1, x0 + 1] * np.ones((wsy, wsx), dtype=int)
        im1 = axs[0].matshow(img1, cmap=cm.gray)
        axs[0].set_title("Original signal")
        for i in range(4):
            xt = x0 + (i % 2) * (wsx + verge) + verge
            yt = y0 + (i // 2) * (wsy + verge) + verge + 30
            axs[0].text(xt, yt, str(enc_orig[y0 + (i // 2), x0 + (i % 2)]), ha="left", va="center", fontsize=6)
        axs[0].xaxis.set_ticks([])
        axs[0].yaxis.set_ticks([])

        x0 = 64 - 1
        y0 = 64 - 1
        img2 = np.abs(coeff)
        img2[y0:(y0 + totaly), x0:(x0 + totalx)] = maxtrafo * np.ones((totaly, totalx), dtype=int)
        img2[(y0 + verge):(y0 + verge + wsy), (x0 + verge):(x0 + verge + wsx)] = coeff[y0, x0] * np.ones((wsy, wsx), dtype=int)
        img2[(y0 + 2 * verge + wsy):(y0 + 2 * verge + 2 * wsy), (x0 + verge):(x0 + verge + wsx)] = coeff[y0 + 1, x0] * np.ones((wsy, wsx), dtype=int)
        img2[(y0 + verge):(y0 + verge + wsy), (x0 + 2 * verge + wsx):(x0 + 2 * verge + 2 * wsx)] = coeff[y0, x0 + 1] * np.ones((wsy, wsx), dtype=int)
        img2[(y0 + 2 * verge + wsy):(y0 + 2 * verge + 2 * wsy), (x0 + 2 * verge + wsx):(x0 + 2 * verge + 2 * wsx)] = coeff[y0 + 1, x0 + 1] * np.ones((wsy, wsx), dtype=int)
        im2 = axs[1].matshow(img2, cmap=cm.gray_r)
        axs[1].set_title("Compressed signal")
        for i in range(4):
            xt = x0 + (i % 2) * (wsx + verge) + verge
            yt = y0 + (i // 2) * (wsy + verge) + verge + 30
            axs[1].text(xt, yt, str(enc_coeff[y0 + (i // 2), x0 + (i % 2)]), ha="left", va="center", fontsize=6)
        axs[1].xaxis.set_ticks([])
        axs[1].yaxis.set_ticks([])

        cbar = plt.colorbar(im1, ax=axs[0], ticks=[0, 256 / 16, 2 * 256 / 16, 255], orientation="vertical")
        cbar.ax.set_yticklabels(["00000000", "00000001", "00000010", "11111111"])
        cbar = plt.colorbar(im2, ax=axs[1], ticks=[0, maxtrafo / 16, 2 * maxtrafo / 16, maxtrafo], orientation="vertical")
        cbar.ax.set_yticklabels(["0", "p10", "p110", "p" + (maxtrafo * "1") + "0"])

        if save_as == "":
            plt.show()
        else:
            assert(save_as.endswith(".svg"))
            print("Processing file '" + save_as + "'...")
            fig.savefig(save_as, format="svg", bbox_inches="tight")

class DemoTrafo:
    @staticmethod
    def visualize(save_as: str = ""):
        #close current opened plot window
        plt.close()
        #generate a simple signal (a combination of two sine waves)
        t = np.linspace(0, 1, 500, endpoint=False)  # Time vector
        signal = 10 * (np.sin(2 * np.pi * 10 * t) + np.sin(2 * np.pi * 50 * t * t))  # Signal

        #perform 1-Level DWT decomposition using Haar wavelet
        cA1, cD1 = pywt.wavedec(signal, 'haar', level=1)
        #perform 3-Level DWT decomposition using Haar wavelet
        #cA3, cD3, cD2, cD1 = pywt.wavedec(signal, 'haar', level=3)

        #quantization with stepsize qs
        qs = 12
        qcA1 = qs * np.round(cA1 / qs)
        qcD1 = qs * np.round(cD1 / qs)

        #perform DWT reconstruction for quantized and non-quantized coeffs
        reco_orig = pywt.waverec([cA1, cD1], wavelet='haar')
        #reco_orig = pywt.waverec([cA3, cD3, cD2, cD1], wavelet='haar')
        reco_quant = pywt.waverec([qcA1, qcD1], wavelet='haar')

        #plot the original signal and the decomposition results
        fig, axs = plt.subplots(5, 1, figsize=(10, 10), layout='constrained')
        c0 = "#006633"
        c1 = "#1f77b4"
        c2 = "#ff7f0e"

        #plot the original signal
        axs[0].plot(t, signal, label="Original", color=c0)
        axs[0].set_title("Original Signal")
        axs[0].set_xlabel("Time (s)")
        axs[0].set_ylabel("Amplitude")

        #plot the approximation coefficients at level 1
        axs[1].plot(cA1, label="non-quantized", color=c1)
        axs[1].plot(qcA1, label="quantized", color=c2)
        axs[1].set_title("DWT: Approximation Coefficients at Level 1")

        #plot the detail coefficients at level 1
        axs[2].plot(cD1, color=c1)
        axs[2].plot(qcD1, color=c2)
        axs[2].set_title("DWT: Detail Coefficients at Level 1")

        #plot signal reconstruction
        axs[3].plot(t, signal, linewidth=1.5, linestyle='dashed', color=c0)
        axs[3].plot(t, reco_orig, linewidth=1.0, color=c1)
        axs[3].set_title("Signal reconstruction")
        axs[3].set_xlabel("Time (s)")
        axs[3].set_ylabel("Amplitude")

        #plot quantized reconstruction
        axs[4].plot(t, signal, linewidth=1.5, linestyle='dashed', color=c0)
        axs[4].plot(t, reco_quant, linewidth=1.0, color=c2)
        axs[4].set_title("Quantized reconstruction")
        axs[4].set_xlabel("Time (s)")
        axs[4].set_ylabel("Amplitude")

        fig.legend(loc='outside right lower')

        if save_as == "":
            plt.show()
        else:
            assert(save_as.endswith(".svg"))
            print("Processing file '" + save_as + "'...")
            fig.savefig(save_as, format="svg", bbox_inches="tight")


class DemoRD:
    @classmethod
    def _collect_data(cls, binImg: str, width: int, height: int, ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list]:
        dist = list()
        bitlen = list()
        qs = list()
        labels = list()
        for predMode in range(0, 5):
            for quantSize in range(8, 65, 8):
                for splitLevel in range(2, 6):
                    dist_val, bitlen_val = run_compression(binImg, width, height, predMode, quantSize, splitLevel)
                    #to make the figure looking nicer, remove higher bitrates
                    if bitlen_val <= 5:
                        dist.append(dist_val)
                        bitlen.append(bitlen_val)
                        qs.append(quantSize)
                        labels.append({"predMode": predMode, "quantSize": quantSize, "splitLevel": splitLevel})
        return np.array(bitlen), np.array(dist), np.array(qs), labels

    @classmethod
    def _get_conv_hull(cls, points: np.ndarray, qs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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

    @classmethod
    def _calc_slopes(cls, points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
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

    @classmethod
    def _interpolate_lambda(cls, vertices: np.ndarray, qs: np.ndarray, slopes: np.ndarray) -> np.ndarray:
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

    @classmethod
    def visualize(cls, binImg: str, width: int, height: int, save_as: str = ""):
        #close current opened plot window
        plt.close()

        bitlen, dist, qs, labels = cls._collect_data(binImg, width, height)
        points = np.column_stack((bitlen, dist))
        vertices, simplices = cls._get_conv_hull(points, qs)
        slopes = cls._calc_slopes(points, vertices)
        lambdas = cls._interpolate_lambda(vertices, qs, slopes)
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

            if save_as == "":
                plt.tight_layout()
                plt.show()
            else:
                fig.set_size_inches(15, 9)
                assert(save_as.endswith(".svg"))
                fname = save_as if figMode == 0 else save_as[:-4] + "_costs.svg"
                print("Processing file '" + fname + "'...")
                fig.savefig(fname, format="svg", bbox_inches="tight", dpi=100)


if __name__ == "__main__":
    img_path = "visuals"
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    img_path += "/"

    svg_steps_exp1 = img_path + "demo_astronaut_pred4_qs8_depth3.svg"
    svg_steps_exp2 = img_path + "demo_astronaut_pred4_qs32_depth3.svg"
    svg_steps_exp3 = img_path + "demo_camera_pred4_qs8_depth2.svg"
    svg_steps_exp4 = img_path + "demo_camera_pred4_qs8_depth5.svg"
    svg_prediction = img_path + "demo_prediction.svg"
    svg_transform  = img_path + "demo_transform.svg"
    svg_encoding   = img_path + "demo_encoding.svg"
    svg_lagrange   = img_path + "demo_RD.svg"
    DemoSteps.visualize("astronaut.bin", 512, 512, predMode=4, quantSize=8, splitLevel=3, save_as=svg_steps_exp1)
    DemoSteps.visualize("astronaut.bin", 512, 512, predMode=4, quantSize=32, splitLevel=3, save_as=svg_steps_exp2)
    DemoSteps.visualize("camera.bin", 512, 512, predMode=4, quantSize=8, splitLevel=2, save_as=svg_steps_exp3)
    DemoSteps.visualize("camera.bin", 512, 512, predMode=4, quantSize=8, splitLevel=5, save_as=svg_steps_exp4)
    DemoPrediction.visualize(save_as=svg_prediction)
    DemoTrafo.visualize(save_as=svg_transform)
    DemoEncoding.visualize("astronaut.bin", 512, 512, predMode=4, quantSize=8, splitLevel=3, save_as=svg_encoding)
    DemoRD.visualize("astronaut.bin", 512, 512, save_as=svg_lagrange)