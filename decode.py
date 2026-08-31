import os
import numpy as np
import struct
from subprocess import run
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

def decode_coef(binstr: list, block_size: int) -> np.ndarray:
    data = np.empty((block_size, block_size))
    i = 0
    for row in range(block_size):
        for col in range(block_size):
            val = 0
            while (binstr[i] != 0):
                val += 1
                i += 1
            i += (1 if val == 0 else 2)
            sgn = 1
            if val > 0:
                sgn = (1 if binstr[i - 1] == 1 else -1)
            data[row, col] = sgn * val
    del binstr[:i]
    return data

def decode_unit(binstr: list, unit_blk_size: int) -> tuple[int, list, list, np.ndarray]:
    dec_depth = binstr[0] * 4 + binstr[1] * 2 + binstr[2]
    del binstr[:3]

    blksize   = unit_blk_size >> dec_depth
    blkcount  = 1 << (dec_depth * 2)
    blkstride = 1 << dec_depth
    dec_pred  = [None] * blkcount
    dec_cut   = [None] * blkcount
    dec_coef  = np.zeros((unit_blk_size, unit_blk_size))

    for bi in range(blkcount):
        dec_pred[bi] = binstr[0] * 2 + binstr[1]
        dec_cut[bi]  = binstr[2] * 2 + binstr[3]
        del binstr[:4]

        row = (bi // blkstride) * blksize
        col = (bi % blkstride) * blksize
        if dec_cut[bi] == 1:
            dec_coef[row:(row + blksize // 2), col:(col + blksize // 2)] = decode_coef(binstr, blksize // 2)
            dec_coef[row:(row + blksize // 2), (col + blksize // 2):(col + blksize)] = decode_coef(binstr, blksize // 2)
            dec_coef[(row + blksize // 2):(row + blksize), col:(col + blksize // 2)] = decode_coef(binstr, blksize // 2)
        elif dec_cut[bi] == 2:
            dec_coef[row:(row + blksize // 2), col:(col + blksize // 2)] = decode_coef(binstr, blksize // 2)
        else:
            dec_coef[row:(row + blksize), col:(col + blksize)] = decode_coef(binstr, blksize)

    return dec_depth, dec_pred, dec_cut, dec_coef

def decode_bitstream(fname: str, unit_blk_size: int) -> tuple[list, list, list, np.ndarray]:
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
    allCoefs  = np.empty((height, width))
    allPreds  = []
    allCuts   = []
    allDepths = []
    for ur in range(height // unit_blk_size):
        for uc in range(width // unit_blk_size):
            currDepth, currPreds, currCuts, currCoef = decode_unit(binstr, unit_blk_size)
            allDepths.append(currDepth)
            allPreds.append(currPreds)
            allCuts.append(currCuts)
            allCoefs[ur * unit_blk_size : (ur + 1) * unit_blk_size, uc * unit_blk_size : (uc + 1) * unit_blk_size] = currCoef
    return allDepths, allPreds, allCuts, allCoefs


if __name__ == "__main__":
    for quantSize in range(4, 25, 4):
        #run encoder and decoder
        os.chdir("comp_demo")
        cmd_comp = "./comp_demo ../astronaut.bin 512 512 " + str(quantSize) + " > /dev/null"
        cmd_reco = "./reco_demo outputs_enc/bitstream.bin > /dev/null"
        run(cmd_comp, shell=True, capture_output=False)
        run(cmd_reco, shell=True, capture_output=False)

        #check if decoder matches encoder
        print("Decoding succesful and valid? ", files_equal("outputs_enc/reco.bin", "outputs_dec/reco.bin"))
        os.chdir("..")

        #check if the Python decoder works correctly
        _, _, _, coef_intern = decode_bitstream("comp_demo/outputs_enc/bitstream.bin", 128)
        bin_file = Path("comp_demo/outputs_enc/coeffs.bin").read_bytes()
        coef_enc = np.reshape(list(struct.unpack('i'*(512 * 512), bin_file)), (512, 512))
        assert((coef_intern == coef_enc).all())