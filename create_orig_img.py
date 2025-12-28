from skimage import data
from skimage.color import rgb2gray
from pathlib import Path
import struct
import numpy as np

#output bitdepth
bitdepth = 8

#create astronaut image from skimage database
grayscale = rgb2gray(data.astronaut())
num_pixel = grayscale.shape[0] * grayscale.shape[1]

print("Creating image 'astronaut.bin' with size " + str(grayscale.shape) + " ...")

grayscale = (np.power(2, bitdepth) - 1) * grayscale.flatten()
grayscale = [np.rint(grayscale[i]).astype(np.int32) for i in range(num_pixel)]

b = struct.pack('i' * num_pixel, *grayscale)
Path("astronaut.bin").write_bytes(b)

#create camera image from skimage database
grayscale = data.camera()
num_pixel = grayscale.shape[0] * grayscale.shape[1]

print("Creating image 'camera.bin' with size " + str(grayscale.shape) + " ...")

grayscale = grayscale.flatten()

b = struct.pack('i' * num_pixel, *grayscale)
Path("camera.bin").write_bytes(b)