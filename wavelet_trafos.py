import numpy as np
import pywt


FilterSet_CDF97 = {
  "h_syn": [ -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 ],
  "g_ana": [ -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 ],
  "h_ana": [  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 ],
  "g_syn": [ -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 ]
}


def convWT(signal: list[float], use_C_implementation=False) -> list:
    scalingFilter  = FilterSet_CDF97["h_ana"]
    hLength        = 9
    waveletFilter  = FilterSet_CDF97["g_ana"]
    gLength        = 7
    n              = len(signal)
    n_out          = n + (0 if use_C_implementation else 8)
    tempbank       = [0.0] * n_out
    filtershift_lo = (((hLength-1) // 2) if use_C_implementation else (hLength - 1))
    filtershift_hi = (((gLength-1) // 2 - 1) if use_C_implementation else (gLength - 1))
    for j in range(0, n_out//2):
        for k in range(0, hLength):
            #symmetric padding before convolution
            idx = 2 * j + k - filtershift_lo
            if idx < 0:
                idx = abs(idx) - 1
            if idx >= n:
                idx = 2 * n - 1 - idx
            tempbank[j] += scalingFilter[k] * signal[idx]
        for k in range(0, gLength):
            idx = 2 * j + k - filtershift_hi
            if idx < 0:
                idx = abs(idx) - 1
            if idx >= n:
                idx = 2 * n - 1 - idx
            tempbank[j + n_out//2] += waveletFilter[k] * signal[idx]
    return tempbank

def invconvWT(trafo: list[float], use_C_implementation=False) -> list:
    scalingFilter  = FilterSet_CDF97["h_syn"]
    hLength        = 7
    waveletFilter  = FilterSet_CDF97["g_syn"]
    gLength        = 9
    n              = len(trafo)
    n_out          = n - (0 if use_C_implementation else 8)
    tempbank       = [0.0] * n_out
    filtershift_lo = (((hLength-1) // 2) if use_C_implementation else -1)
    filtershift_hi = (((gLength-1) // 2) if use_C_implementation else 0)
    for j in range(0, n_out):
        for k in range(0, hLength):
            idx = j + k - filtershift_lo
            if idx % 2 == 0:
                if idx < 0:
                    idx = abs(idx) - 1
                if idx >= n:
                    idx = 2 * n - 1 - idx
                tempbank[j] += scalingFilter[k] * trafo[(idx//2)]
        for k in range(0, gLength):
            idx = j + k - filtershift_hi
            if idx % 2 == 1:
                if idx < 0:
                    idx = abs(idx) - 1
                if idx >= n:
                    idx = 2 * n - 1 - idx
                tempbank[j] += waveletFilter[k] * trafo[(idx//2 + n//2)]
    return tempbank

def generate_test_signal() -> tuple[np.ndarray, np.ndarray]:
    t = np.linspace(0, 1, 500, endpoint=False)  # Time vector
    f = 10 * (np.sin(2 * np.pi * 10 * t) + np.sin(2 * np.pi * 50 * t * t)) # Signal
    return t, f


if __name__ == "__main__":
    x, y = generate_test_signal()
    assert(len(y) >= 16)

    cA1, cD1 = pywt.wavedec(y, "bior4.4", level=1, mode="symmetric")
    reco_signal = pywt.waverec([cA1, cD1], wavelet="bior4.4", mode="symmetric")

    print("Filter set:")
    print(pywt.Wavelet("bior4.4").filter_bank)
    print()

    coefs = convWT(y)

    coef_low = coefs[:(len(coefs)//2)]
    coef_high = coefs[(len(coefs)//2):]
    assert(len(cA1)==len(coef_low))
    assert(np.isclose(coef_low, cA1, atol=0.0001).all())
    assert(np.isclose(coef_high, cD1, atol=0.0001).all())

    rec = invconvWT(coefs)

    assert(len(reco_signal) == len(rec))
    assert(np.isclose(reco_signal, rec, atol=0.0001).all())

    coefs_c = convWT(y, use_C_implementation=True)
    rec_c = invconvWT(coefs_c, use_C_implementation=True)

    assert(np.isclose(reco_signal[5:-5], rec_c[5:-5], atol=0.0005).all())
    mse = np.square(np.subtract(rec_c, reco_signal)).mean()
    print("MSE: " + str(mse) + "\n")