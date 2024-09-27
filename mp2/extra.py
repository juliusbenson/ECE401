'''
If you finish this module, you can submit it for extra credit.
'''
import numpy as np
import submitted
import cmath

def rate_conversion(highrate_signal, highrate, lowrate):
    '''
    lowrate_signal, coefficients = rate_conversion(highrate_signal, highrate, lowrate)
    
    highrate_signal (real array, length N_0*T) - one period of the high-rate signal
    highrate (scalar) - the high sampling rate, expressed in samples/second
    lowrate (scalar) - the low sampling rate, expressed in samples/second
    lowrate_signal (real array, length N_0) - one period of the low-rate signal
    coefficients (complex array) - Fourier coefficients X_{-int((N0-1)/2)} through X_{int((N0-1)/2)}.
    
    This should perform a Fourier series analysis of the high-rate signal.
    Then resynthesize the low-rate signal using Fourier synthesis, with only the terms
    that are still valid at the low-rate, i.e., X_{-int((N0-1)/2)} through X_{int((N0-1)/2)}.

    You'll want to use a complex-valued lowrate_signal during Fourier synthesis,
    then return its real part.
    '''

    N_0T = len(highrate_signal)
    print(f'N_0T = {N_0T}')
    T = highrate / lowrate
    print(f'T = {T}')
    N_0 = int(np.floor(N_0T / T))
    print(f'N_0 = {N_0}')
    coefficients = submitted.fourier_analysis(highrate_signal,int(np.ceil((N_0-1)/2)))
    print(f'{len(coefficients)} coefficients: {coefficients}')
    lowrate_signal = np.zeros(int(N_0))
    for n,losample in enumerate(lowrate_signal):
        lowrate_signal[n] = np.sum([(2*X_k*cmath.rect(1,2*cmath.pi*k*n/N_0)) for k,X_k in enumerate(coefficients)]).real

    coefficients = np.array(list(reversed(coefficients)) + coefficients[1:])
    print(f'{len(coefficients)} coefficients returned: {coefficients}')

    return lowrate_signal, coefficients
    raise RuntimeError("You need to write this!")
