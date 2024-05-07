from skimage import data
from skimage.color import rgb2gray
from pathlib import Path
import struct
import numpy as np

##############################################
#create grayscale image from skimage database#
##############################################
bitdepth = 8

original = data.astronaut()
grayscale = rgb2gray(original)
num_pixel = grayscale.shape[0] * grayscale.shape[0]
print("width x height of picture: " + str(grayscale.shape))

grayscale = (np.power(2, bitdepth) - 1) * grayscale.flatten()
grayscale = [np.rint(grayscale[i]).astype(np.int32) for i in range(num_pixel)]
#print(grayscale)

b = struct.pack('i'*num_pixel, *grayscale)
Path('astronaut.bin').write_bytes(b)


#########################################
#create mathematical grayscale image 8x8#
#########################################
grayscale2 = []
for i in range(8):
    for j in range(8):
      grayscale2.append(50+i+j+4*(3-i)*(5-j)+j*j)
#print(grayscale2)

b = struct.pack('i'*64, *grayscale2)
Path('art8x8.bin').write_bytes(b)