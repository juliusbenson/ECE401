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
    T = highrate / lowrate
    N_0 = N_0T / T
    coefficients = submitted.fourier_analysis(highrate_signal,int((N_0-1)//2)+1)
    print(f'{len(coefficients)} coefficients: {coefficients}')
    lowrate_signal = np.zeros(int(N_0))
    for n,losample in enumerate(lowrate_signal):
        # toSum = 0 + 0j
        # for k,X_k in enumerate(coefficients): # these will only be the positive coefficients, but we know that the negative coefficients have the same frequency (but negative), and the same phase and amplitude
        # # for k,X_k in enumerate(coefficients):
        #     # X_k = coefficients[k]
        #     print(f'coefficient X_{k} = {X_k}')
        #     toSum += cmath.rect(X_k,2*cmath.pi*k*n/N_0) + cmath.rect(X_k,2*cmath.pi*(-k)*n/N_0) # doing two symmetrical coefficients at once, we're all adding them together anyway
        # if toSum.imag != 0: raise Exception("fourier sum contains imaginary part") # if we get here, something's wrong in the code. we should never hit this.
        # lowrate_signal[n] = toSum.real
        lowrate_signal[n] = np.sum([cmath.rect(X_k,2*cmath.pi*k*n/N_0).real for k,X_k in enumerate(coefficients)])

    coefficients = np.array(list(reversed(coefficients)) + coefficients[1:])
    print(len(coefficients))
    print(coefficients)

    return lowrate_signal, coefficients
    raise RuntimeError("You need to write this!")
