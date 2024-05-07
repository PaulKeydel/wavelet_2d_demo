import sys
from pathlib import Path
import struct
import numpy as np
import matplotlib.pyplot as plt
from pylab import cm
from matplotlib.colors import SymLogNorm

assert(len(sys.argv) == 3)
width = int(sys.argv[1])
height = int(sys.argv[2])
num_pixel = width * height

orig_bytes = Path('orig.bin').read_bytes()
orig_data = struct.unpack('i'*num_pixel, orig_bytes)
orig_data = np.reshape(orig_data, (width,height))

reco_bytes = Path('reco.bin').read_bytes()
reco_data = struct.unpack('i'*num_pixel, reco_bytes)
reco_data = np.reshape(reco_data, (width,height))

coeff_bytes = Path('coeffs.bin').read_bytes()
coeff_data = struct.unpack('i'*num_pixel, coeff_bytes)
coeff_data = np.reshape(coeff_data, (width,height))
coeff_data = abs(coeff_data)

coeffq_bytes = Path('coeffs_quant.bin').read_bytes()
coeffq_data = struct.unpack('i'*num_pixel, coeffq_bytes)
coeffq_data = np.reshape(coeffq_data, (width,height))
coeffq_data = abs(coeffq_data)

fig, axs = plt.subplots(2, 2)
axs[0, 0].matshow(orig_data)
axs[0, 0].set_title('orig')

pcoeff = axs[0, 1].matshow(coeff_data, cmap=cm.gray_r, norm=SymLogNorm(linthresh=4, base=2))
fig.colorbar(pcoeff, ax=axs[0, 1])
axs[0, 1].set_title('abs(coeff)')

axs[1, 0].matshow(reco_data)
axs[1, 0].set_title('reco')

pcoeffq = axs[1, 1].matshow(coeffq_data, cmap=cm.gray_r, norm=SymLogNorm(linthresh=4, base=2))
fig.colorbar(pcoeffq, ax=axs[1, 1])
axs[1, 1].set_title('quant abs(coeff)')
#pdiff = axs[1, 1].matshow(reco_data - orig_data, cmap=cm.gray_r)
#fig.colorbar(pdiff, ax=axs[1, 1])
#axs[1, 1].set_title('reco - orig')

plt.show()