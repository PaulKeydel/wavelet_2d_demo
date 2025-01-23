from skimage import data
from skimage.color import rgb2gray
from pathlib import Path
import struct
import numpy as np

#output bitdepth
bitdepth = 8

##############################################
#create grayscale image from skimage database#
##############################################
grayscale = rgb2gray(data.astronaut())
num_pixel = grayscale.shape[0] * grayscale.shape[0]

print("Creating image 'astronaut.bin' with size " + str(grayscale.shape) + " ...")

grayscale = (np.power(2, bitdepth) - 1) * grayscale.flatten()
grayscale = [np.rint(grayscale[i]).astype(np.int32) for i in range(num_pixel)]
#print(grayscale)

b = struct.pack('i' * num_pixel, *grayscale)
Path('astronaut.bin').write_bytes(b)


#########################################
#create mathematical grayscale image 1024x1024
#########################################
grayscale2 = []
num_pixel = 1024 * 1024

print("Creating image 'radial1024.bin' with size (1024, 1024) ...")

for i in range(1024):
  for j in range(1024):
    pixel = np.rint(np.sqrt(i ** 2 + j ** 2) * (2 ** bitdepth - 1) / 1023 / np.sqrt(2)).astype('int')
    grayscale2.append(pixel)

b = struct.pack('i' * num_pixel, *grayscale2)
Path('radial1024.bin').write_bytes(b)