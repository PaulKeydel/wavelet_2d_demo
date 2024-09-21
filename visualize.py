import sys
from pathlib import Path
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import cm
from matplotlib.colors import SymLogNorm

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

assert(len(sys.argv) == 3)
width = int(sys.argv[1])
height = int(sys.argv[2])
num_pixel = width * height

orig_bytes = Path('orig.bin').read_bytes()
orig_data = struct.unpack('i'*num_pixel, orig_bytes)
entr_orig = entropy(orig_data)
orig_data = np.reshape(orig_data, (width,height))

reco_bytes = Path('reco.bin').read_bytes()
reco_data = struct.unpack('i'*num_pixel, reco_bytes)
reco_data = np.reshape(reco_data, (width,height))

coeff_bytes = Path('coeffs.bin').read_bytes()
coeff_data = struct.unpack('i'*num_pixel, coeff_bytes)
entr_coeff = entropy(coeff_data)
coeff_data = np.reshape(coeff_data, (width,height))
coeff_data = abs(coeff_data)

resi_bytes = Path('resi.bin').read_bytes()
resi_data = struct.unpack('i'*num_pixel, resi_bytes)
resi_data = np.reshape(resi_data, (width,height))
resi_data = abs(resi_data)

fig, axs = plt.subplots(2, 2)
axs[0, 0].matshow(orig_data)
axs[0, 0].set_title('orig')

axs[1, 0].matshow(reco_data)
axs[1, 0].set_title('reco')

presi = axs[0, 1].matshow(resi_data, cmap=cm.gray_r, norm=SymLogNorm(linthresh=4, base=2))
fig.colorbar(presi, ax=axs[0, 1])
axs[0, 1].set_title('abs(resi)')

pcoeff = axs[1, 1].matshow(coeff_data, cmap=cm.gray_r, norm=SymLogNorm(linthresh=2, base=2))
fig.colorbar(pcoeff, ax=axs[1, 1])
axs[1, 1].set_title('abs(quant coeff)')

txt = "Entropy original image: " + "{:.3f}".format(entr_orig) + "   /   entropy transformed image: " + "{:.3f}".format(entr_coeff)
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

plt.show()