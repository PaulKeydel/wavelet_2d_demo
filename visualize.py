import sys
from pathlib import Path
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pylab import cm

def collect_data(num_pixel: int) -> tuple[list, list, list, list]:
    orig_bytes = Path('orig.bin').read_bytes()
    orig_data = list(struct.unpack('i'*num_pixel, orig_bytes))

    reco_bytes = Path('reco.bin').read_bytes()
    reco_data = list(struct.unpack('i'*num_pixel, reco_bytes))

    coeff_bytes = Path('coeffs.bin').read_bytes()
    coeff_data = list(struct.unpack('i'*num_pixel, coeff_bytes))

    resi_bytes = Path('resi.bin').read_bytes()
    resi_data = list(struct.unpack('i'*num_pixel, resi_bytes))

    return orig_data, reco_data, resi_data, coeff_data

def entropy(message):
    n_labels = len(message)
    if n_labels <= 1:
        return 0
    value, counts = np.unique(message, return_counts=True)
    probs = counts / n_labels
    n_classes = len(value)
    if n_classes <= 1:
        return 0
    ent = 0.
    for i in probs:
        ent -= i * math.log(i, 2)
    return ent

def make_2d_image(src: list, width: int, height: int, flip_scale = False) -> np.ndarray:
    ptp = np.ptp(src)
    scaled_8bit = []
    if not flip_scale:
        min_val = min(src)
        scaled_8bit = [255 * (src[i] - min_val) / ptp for i in range(len(src))]
    else:
        max_val = max(src)
        scaled_8bit = [255 * (max_val - src[i]) / ptp for i in range(len(src))]
    return np.reshape(scaled_8bit, (width, height))

def plot_bins(width: int, height: int, save_as: str | None):
    orig_data, reco_data, resi_data, coeff_data = collect_data(width * height)

    print("  orig: max=" + str(max(orig_data)) + " min=" + str(min(orig_data)))
    print("  reco: max=" + str(max(reco_data)) + " min=" + str(min(reco_data)))
    print("  coef: max=" + str(max(coeff_data)) + " min=" + str(min(coeff_data)))
    print("  resi: max=" + str(max(resi_data)) + " min=" + str(min(resi_data)))

    entr_orig = entropy(orig_data)
    entr_coeff = entropy(coeff_data)

    coeff_data = list(map(abs, coeff_data))
    resi_data = list(map(abs, resi_data))

    clip_x0 = 0
    clip_x1 = 256
    clip_y0 = 0
    clip_y1 = 256
    img1 = make_2d_image(orig_data, width, height)
    img2 = make_2d_image(reco_data, width, height)
    img3 = make_2d_image(resi_data, width, height, flip_scale=True)[clip_y0:clip_y1, clip_x0:clip_x1]
    img4 = make_2d_image(coeff_data, width, height, flip_scale=True)[clip_y0:((clip_y0 + clip_y1) // 2), clip_x0:((clip_x0 + clip_x1) // 2)]

    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    axs[0, 0].matshow(img1, cmap=cm.gray)
    axs[0, 0].set_title("Original")
    axs[0, 0].xaxis.set_ticks_position("bottom")
    if clip_x0 != 0 or clip_x1 != width or clip_y0 != 0 or clip_y1 != height:
        rect = patches.Rectangle((clip_x0, clip_y0), clip_x1 - clip_x0, clip_y1 - clip_y0, linewidth=1, edgecolor='r', facecolor='none')
        axs[0, 0].add_patch(rect)

    axs[1, 0].matshow(img2, cmap=cm.gray)
    axs[1, 0].set_title("Reconstructed")
    axs[1, 0].xaxis.set_ticks([])
    axs[1, 0].yaxis.set_ticks([])

    axs[0, 1].matshow(img3, cmap=cm.gray)
    axs[0, 1].set_title("Residual (absolute values)")
    axs[0, 1].xaxis.set_ticks([])
    axs[0, 1].yaxis.set_ticks([])
    if clip_x0 != 0 or clip_x1 != width or clip_y0 != 0 or clip_y1 != height:
        rect = patches.Rectangle((clip_x0, clip_y0), (clip_x1 - clip_x0) // 2, (clip_y1 - clip_y0) // 2, linewidth=1, edgecolor='r', facecolor='none')
        axs[0, 1].add_patch(rect)

    axs[1, 1].matshow(img4, cmap=cm.gray)
    axs[1, 1].set_title("Transformed and quantized")
    axs[1, 1].xaxis.set_ticks([])
    axs[1, 1].yaxis.set_ticks([])

    txt = "Entropy original image: " + "{:.3f}".format(entr_orig) + "   /   entropy transformed image: " + "{:.3f}".format(entr_coeff)
    plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

    if save_as is None:
        plt.show()
    else:
        assert(save_as.endswith(".svg"))
        fig.savefig(save_as, format="svg", bbox_inches="tight")

if __name__ == "__main__":
    assert(len(sys.argv) == 3 or len(sys.argv) == 4)
    width = int(sys.argv[1])
    height = int(sys.argv[2])

    if (len(sys.argv) == 3):
        plot_bins(width, height, save_as=None)
    else:
        plot_bins(width, height, save_as=sys.argv[3])