import os
import numpy as np
import struct
from pathlib import Path

def files_equal(path1, path2, chunk=8192): 
    with open(path1, 'rb') as f1, open(path2, 'rb') as f2: 
        while True: 
            b1 = f1.read(chunk) 
            b2 = f2.read(chunk) 
            if b1 != b2: 
                return False 
            if not b1: #both reached EOF 
                return True

def decode_unit(binstr: list, unit_blk_size: int) -> tuple[int, int, np.ndarray, int]:
    dec_depth = binstr[0] * 4 + binstr[1] * 2 + binstr[2]
    dec_pred  = binstr[3] * 2 + binstr[4]
    dec_coef  = np.empty((unit_blk_size, unit_blk_size))
    i         = 5
    i_last    = 5
    for row in range(unit_blk_size):
        for col in range(unit_blk_size):
            while (binstr[i] != 0):
                i += 1
            dc = i - i_last
            if dc > 0:
                i += 1
                dc *= (1 if binstr[i] == 1 else -1)
            i += 1
            dec_coef[row, col] = dc
            i_last = i
    if i % 8 != 0:
        i += 8 - (i % 8)
    assert(i % 8 == 0)
    return dec_depth, dec_pred, dec_coef, i

def decode_full(fname: str, unit_blk_size: int) -> tuple[list, list, np.ndarray]:
    binstr = []
    with open(fname, "rb") as f:
        coding_params = f.read(4) # 4 header bytes for width, height and quant step-size
        width = int.from_bytes(coding_params[:3], "big") >> 12
        height = int.from_bytes(coding_params[:3], "big") & 4095
        byte = f.read(1)
        while byte:
            byteval = int.from_bytes(byte, "big")
            for i in range(8):
                binstr.append((byteval >> (7 - i)) & 1)
            byte = f.read(1)
    #decode unit-wise
    data_coef  = np.empty((height, width))
    data_pred  = []
    data_depth = []
    offset = 0
    for ur in range(height // unit_blk_size):
        for uc in range(width // unit_blk_size):
            depth_read, pred_read, coef_read, bits_read = decode_unit(binstr[offset:], unit_blk_size)

            data_depth.append(depth_read)
            data_pred.append(pred_read)

            data_coef[ur * unit_blk_size : (ur + 1) * unit_blk_size, uc * unit_blk_size : (uc + 1) * unit_blk_size] = coef_read

            offset += bits_read
    return data_depth, data_pred, data_coef


if __name__ == "__main__":
    for quantSize in range(4, 25, 4):
        #run encoder and decoder
        os.chdir("comp_demo")
        os.system("./comp_demo ../astronaut.bin 512 512 " + str(quantSize) + " > /dev/null")
        os.system("./reco_demo outputs_enc/bitstream.bin > /dev/null")

        #check if decoder matches encoder
        print("Decoding succesful and valid? ", files_equal("outputs_enc/reco.bin", "outputs_dec/reco.bin"))
        os.chdir("..")

        #check if the Python decoder works correctly
        _, _, coef_intern = decode_full("comp_demo/outputs_enc/bitstream.bin", 128)
        buffer = Path("comp_demo/outputs_enc/coeffs.bin").read_bytes()
        coef_enc = np.reshape(list(struct.unpack('i'*(512 * 512), buffer)), (512, 512))
        assert((coef_intern == coef_enc).all())