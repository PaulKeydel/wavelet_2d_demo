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

def decode_unit(binstr: list, unit_blk_size: int) -> tuple[int, list, list, np.ndarray, int]:
    dec_depth = binstr[0] * 4 + binstr[1] * 2 + binstr[2]
    blksize   = unit_blk_size >> dec_depth
    blkcount  = 1 << (dec_depth * 2)
    blkstride = 1 << dec_depth
    dec_pred  = [None] * blkcount
    dec_cut   = [None] * blkcount
    dec_coef  = np.empty((unit_blk_size, unit_blk_size))
    i         = 3
    i_last    = 3
    for bi in range(blkcount):
        dec_pred[bi] = binstr[i + 0] * 2 + binstr[i + 1]
        dec_cut[bi]  = binstr[i + 2] * 2 + binstr[i + 3]
        i += 4
        i_last += 4
        blkrow = (bi // blkstride) * blksize
        blkcol = (bi % blkstride) * blksize
        for row in range(blkrow, blkrow + blksize):
            for col in range(blkcol, blkcol + blksize):
                if ((dec_cut[bi] == 1 and (row % blksize >= blksize // 2 and col % blksize >= blksize // 2)) or
                    (dec_cut[bi] == 2 and (row % blksize >= blksize // 2 or col % blksize >= blksize // 2))):
                    dec_coef[row, col] = 0
                    continue
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
    return dec_depth, dec_pred, dec_cut, dec_coef, i

def decode_full(fname: str, unit_blk_size: int) -> tuple[list, list, list, np.ndarray]:
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
    offset = 0
    for ur in range(height // unit_blk_size):
        for uc in range(width // unit_blk_size):
            currDepth, currPreds, currCuts, currCoef, bits_read = decode_unit(binstr[offset:], unit_blk_size)

            allDepths.append(currDepth)
            allPreds.append(currPreds)
            allCuts.append(currCuts)
            allCoefs[ur * unit_blk_size : (ur + 1) * unit_blk_size, uc * unit_blk_size : (uc + 1) * unit_blk_size] = currCoef

            offset += bits_read
    return allDepths, allPreds, allCuts, allCoefs

def write_cfg_overview(fname: str):
    bestDepths, bestPreds, bestCuts, coef_intern = decode_full("comp_demo/outputs_enc/bitstream.bin", 128)
    #check if the Python decoder works correctly
    buffer = Path("comp_demo/outputs_enc/coeffs.bin").read_bytes()
    coef_enc = np.reshape(list(struct.unpack('i'*(512 * 512), buffer)), (512, 512))
    assert((coef_intern == coef_enc).all())
    #write unit-related config overview in file
    with open(fname, "w") as f:
        f.write("Unit-related coding parameters\n")
        f.write("------------------------------\n")
        for i in range(len(bestDepths)):
            f.write("depth: " + str(bestDepths[i]) + "\n")
            f.write("  preds: " + str(bestPreds[i]) + "\n")
            f.write("  cuts:  " + str(bestCuts[i]) + "\n")
        f.close()


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

    write_cfg_overview("decoded_config.txt")