import sys
from pathlib import Path
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import cm

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

def normalize_to_8bit(vec, reverse = False):
    ptp_ = np.ptp(vec)
    if (not reverse):
        min_ = min(vec)
        return [255 * (vec[i] - min_) / ptp_ for i in range(len(vec))]
    else:
        max_ = max(vec)
        return [255 * (max_ - vec[i]) / ptp_ for i in range(len(vec))]

assert(len(sys.argv) == 3 or len(sys.argv) == 4)
width = int(sys.argv[1])
height = int(sys.argv[2])
num_pixel = width * height

orig_bytes = Path('orig.bin').read_bytes()
orig_data = struct.unpack('i'*num_pixel, orig_bytes)
print("  orig: max=" + str(max(orig_data)) + " min=" + str(min(orig_data)))
entr_orig = entropy(orig_data)
orig_data = normalize_to_8bit(orig_data)
orig_data = np.reshape(orig_data, (width,height))

reco_bytes = Path('reco.bin').read_bytes()
reco_data = struct.unpack('i'*num_pixel, reco_bytes)
print("  reco: max=" + str(max(reco_data)) + " min=" + str(min(reco_data)))
reco_data = normalize_to_8bit(reco_data)
reco_data = np.reshape(reco_data, (width,height))

coeff_bytes = Path('coeffs.bin').read_bytes()
coeff_data = struct.unpack('i'*num_pixel, coeff_bytes)
print("  coeff: max=" + str(max(coeff_data)) + " min=" + str(min(coeff_data)))
assert(min(coeff_data) == 0)
entr_coeff = entropy(coeff_data)
coeff_data = normalize_to_8bit(coeff_data, reverse=True)
coeff_data = np.reshape(coeff_data, (width,height))

resi_bytes = Path('resi.bin').read_bytes()
resi_data = struct.unpack('i'*num_pixel, resi_bytes)
print("  resi: max=" + str(max(resi_data)) + " min=" + str(min(resi_data)))
resi_data = normalize_to_8bit(list(map(abs, resi_data)), reverse=True)
#print("  resi: max=" + str(max(resi_data)) + " min=" + str(min(resi_data)))
resi_data = np.reshape(resi_data, (width,height))

fig, axs = plt.subplots(2, 2, figsize=(8, 8))

axs[0, 0].matshow(orig_data, cmap=cm.gray)
axs[0, 0].set_title('Original')
axs[0, 0].xaxis.set_ticks([])

axs[1, 0].matshow(reco_data, cmap=cm.gray)
axs[1, 0].set_title('Reconstructed')
axs[1, 0].xaxis.set_ticks_position('bottom')

axs[0, 1].matshow(resi_data, cmap=cm.gray)
axs[0, 1].set_title('Residual (absolute values)')
axs[0, 1].xaxis.set_ticks([])
axs[0, 1].yaxis.set_ticks([])

axs[1, 1].matshow(coeff_data, cmap=cm.gray)
axs[1, 1].xaxis.set_ticks_position('bottom')
axs[1, 1].yaxis.set_ticks([])
axs[1, 1].set_title('Transformed and quantized')

txt = "Entropy original image: " + "{:.3f}".format(entr_orig) + "   /   entropy transformed image: " + "{:.3f}".format(entr_coeff)
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

if (len(sys.argv) == 3):
    plt.show()
else:
    plt.savefig(sys.argv[3], format="svg", bbox_inches="tight")