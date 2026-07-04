import os
import numpy as np
import struct
from pathlib import Path

def decode_unit(binstr: list, nitems: int):
    dec_depth = binstr[0] * 4 + binstr[1] * 2 + binstr[2]
    dec_pred  = binstr[3] * 2 + binstr[4]
    dec_coef  = []
    i = 5
    count_ones = 0
    while len(dec_coef) < nitems:
        if binstr[i] == 0:
            dc = count_ones
            if count_ones > 0:
                i += 1
                dc *= (1 if binstr[i] == 1 else -1)
            dec_coef.append(dc)
            count_ones = 0
        else:
            count_ones += 1
        i += 1
    assert(len(dec_coef) == nitems)
    if i % 8 != 0:
        i += 8 - (i % 8)
    assert(i % 8 == 0)
    return dec_depth, dec_pred, dec_coef, i

def decode_full(fname: str, width: int, height: int, unit_blk_size: int) -> tuple[list, list, np.ndarray]:
    binstr = []
    with open(fname, "rb") as f:
        byte = f.read(1)
        while byte:
            byteval = int.from_bytes(byte, "big")
            for i in range(8):
                binstr.append((byteval >> (7 - i)) & 1)
            byte = f.read(1)
    #decode binary stream now
    full_coef  = np.empty((height, width))
    full_pred  = []
    full_depth = []
    offset = 0
    nitems = unit_blk_size * unit_blk_size
    for ur in range(height // unit_blk_size):
        for uc in range(width // unit_blk_size):
            depth_read, pred_read, coeff_read, bits_read = decode_unit(binstr[offset:], nitems)
            full_depth.append(depth_read)
            full_pred.append(pred_read)
            full_coef[ur * unit_blk_size : (ur + 1) * unit_blk_size, uc * unit_blk_size : (uc + 1) * unit_blk_size] = np.reshape(coeff_read, (unit_blk_size, unit_blk_size))
            offset += bits_read
    return full_depth, full_pred, full_coef


if __name__ == "__main__":
    for quantSize in range(4, 25, 4):
        os.system("./comp_demo astronaut.bin 512 512 " + str(quantSize) + " > /dev/null")

        _, _, full_coef = decode_full("bitstream.bin", 512, 512, 128)

        buffer = Path("coeffs.bin").read_bytes()
        coeff = np.reshape(list(struct.unpack('i'*(512 * 512), buffer)), (512, 512))

        print("Decoding succesful and valid?", (full_coef == coeff).all())