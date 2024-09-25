'''
If you finish this module, you can submit it for extra credit.
'''
import numpy as np

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
    raise RuntimeError("You need to write this!")
