'''
If you finish this module, you can submit it for extra credit.
'''
import numpy as np
import cmath

def fourier_series_coefficients(x, F0, Fs, N):
    '''
    phasors = fourier_series_coefficients(x, F0, Fs, N)
    
    x (numpy array) - the input signal
    F0 (scalar) - the fundamental frequency, in Hertz
    Fs (scalar) - sampling frequency, in samples/second
    N (scalar) - number of harmonics to measure
    phasors (length-N array) - Fourier series coefficients
    
    This should compute only the positive-frequency Fourier series coefficients
    (for k=1 through k=N).  Assume that the first sample of the input signal, x[0],
    is at time 0; the second sample, x[1], is at time 1/Fs, and so on.
    
    Instead of averaging over one period, you should average over the 
    whole signal length.  So you should multiply x by the complex 
    conjugate of each harmonic, then add over the whole signal length,
    then divide by the length of the signal.
    '''

    xks = []
    for k in range(N+1):
        xns = []
        for n,xn in enumerate(x):
            tosum = cmath.rect(xn,-2 * np.pi * k * n * F0 / Fs)
            xns.append(tosum)
        xk = np.sum(xns) / len(x)
        print(xk)
        xks.append(xk)
    return xks[1:]
    # raise RuntimeError("You need to write this!")
