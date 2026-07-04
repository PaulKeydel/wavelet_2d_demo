import numpy as np
from visualize import (run_compression, load_binaries)

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


for quantSize in range(4, 25, 4):
    run_compression("astronaut.bin", 512, 512, quantSize)

    binstr = []
    with open("bitstream.bin", "rb") as f:
        byte = f.read(1)
        while byte:
            byteval = int.from_bytes(byte, "big")
            for i in range(8):
                binstr.append((byteval >> (7 - i)) & 1)
            byte = f.read(1)

    full_coef = np.empty((512, 512))
    offset = 0
    for ur in range(4):
        for uc in range(4):
            depth_read, pred_read, coeff_read, bits_read = decode_unit(binstr[offset:], 128 * 128)
            #print("depth: " + str(depth_read))
            #print("pred: " + str(pred_read))
            full_coef[ur * 128 : (ur + 1) * 128, uc * 128 : (uc + 1) * 128] = np.reshape(coeff_read, (128, 128))
            offset += bits_read

    _, _, _, _, coeff = load_binaries(512, 512)
    print("Decoding succesful and valid?", (full_coef == coeff).all())