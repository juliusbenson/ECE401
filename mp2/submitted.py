'''
This is the module you'll submit to the autograder.

There are several function definitions, here, that raise RuntimeErrors.  You should replace
each "raise RuntimeError" line with a line that performs the function specified in the
function's docstring.
'''

import numpy as np
import cmath

def sinusoid(phasor:complex, frequency:float, duration:float, samplerate:float) -> tuple[list[float], list[float]]:
    '''
    timeaxis, signal = sinusoid(phasor, frequency, duration, samplerate)
    Generate a sinusoid.

    phasor (complex scalar) - magnitude times e^{j phase}
    frequency (real scalar) - frequency of the sinusoid, in Hertz
    duration (real scalar) - duration, in seconds
    samplerate (real scalar) - sampling rate, in samples/second
    timeaxis (array) - sample times, from 0 to duration, including endpoints
    signal (array) - the generated sinusoid, length = int(duration*samplerate+1)
    '''

    samples = (int)(samplerate*duration) # ( samples / second ) * seconds = samples
    timeaxis = [n/samplerate for n in range(0,samples+1)] # inluding endpoint
    signal = [(phasor * cmath.rect(1,2*cmath.pi*t*frequency)).real for t in timeaxis]
    return timeaxis, signal
    
def compute_aliasing(frequencies:list[float], phasors:list[complex], samplerates:list[float]) -> tuple[list[float], list[complex]]:
    '''
    aliased_freqs, aliased_phasors = compute_aliasing(frequencies, phasors, samplerates)
    Find the frequency and phasor of sinusoid aliases.  All arguments should have same length.

    frequencies (real array) - frequencies of the sinusoids, in Hertz
    phasors (complex array) - magnitudes times e^{j phases}
    samplerates (real array) - sampling rates, in samples/second
    aliased_freqs (real array)  - frequencies at which sinusoids seems to occur, in Hertz
    aliased_phasors (complex array) - phasors with which sinusoids seems to occur
    '''

    aliased_phasors = []
    for k,_ in enumerate(frequencies):
        f =  frequencies[k]
        z =      phasors[k]
        Fs = samplerates[k]

        if f / Fs < round (f / Fs):
            za = z.conjugate()
        else:
            za = z

        aliased_phasors.append(za)

    aliased_freqs = [min(frequencies[k] % samplerates[k], (samplerates[k]-frequencies[k]) % samplerates[k]) for k in range(len(frequencies))]

    return aliased_freqs, aliased_phasors

def fourier_analysis(signal : list[float], number_of_coefficients : int) -> list[complex]:
    '''
    coefficients = fourier_analysis(signal, number_of_coefficients)
    Find the Fourier series coefficients using the discrete-time Fourier analysis formula.

    signal (array of length N_0) = one period of the signal
    number_of_coefficients (scalar) = number of coefficients to compute, starting with X_0
    coefficients (array of length=number_of_coefficients) = X_0 through X_{number_of_coefficients-1}
    '''

    N_0 = len(signal)
    coefficients = []

    for k in range(number_of_coefficients):
        X_k = np.sum([(cmath.rect(signal[n],-2*cmath.pi*k*n/N_0)) for n in range(N_0)]) / N_0
        coefficients.append(X_k)

    return coefficients

def interpolate(lowrate_signal : list[float], T : float, kernel_timeaxis : list[float], kernel : list[float]) -> list[float]:
    '''
    highrate_signal = interpolate(lowrate_signal, T, kernel_timeaxis, kernel)
    Use lowrate-to-highrate conversion to simulate discrete-to-continuous conversion.

    lowrate_signal (length-N array) - the lowrate signal
    T (scalar) - ratio of highrate/lowrate, i.e., number of output samples per input sample
    kernel_timeaxis (array) - sample times of the kernel, at the highrate
    kernel (array) - the interpolation kernel.  length(kernel)==length(kernel_timeaxis).
    highrate_signal (length-N*T array) - the highrate signal
    '''
    def h(time) -> float: # address kernel values by time, not by array index
        kernel[kernel_timeaxis.index(time)]

    highrate_signal = np.zeros(T*len(lowrate_signal))

    for t,hisample in enumerate(highrate_signal):
        for n,losample in enumerate(lowrate_signal):
            tosum = losample * h(t - (n*T))
        hisample = np.sum(tosum)

    return highrate_signal
    raise RuntimeError("You need to write this part!")

def rectangle(T : float) -> tuple[list[int], list[float]]:
    '''
    timeaxis, h = rectangle(T)
    Return a rectangle function of length T.

    T (scalar) - length, in samples
    timeaxis (length-T array) - sample indices, from 0 to T-1, corresponding to h
    h (length-T array) - the rectangle function
    '''
    raise RuntimeError("You need to write this part!")

def triangle(T : int) -> tuple[list[int], list[float]]:
    '''
    timeaxis, h = triangle(T)
    Return a triangle function of length 2*T-1.

    T (scalar) - length of each side of the triangle, in samples
    timeaxis (array, length 2*T-1) - sample indices, from -(T-1) through (T-1)
    h (array, length 2*T-1) - the triangle function, 1 - abs(timeaxis)/T
    '''
    timeaxis = (list)(range(-T+1,T))
    h = []

    for n in timeaxis:
        hn = 1 - np.abs(n) / T
        h.append(hn)

    print(timeaxis)
    return timeaxis, h
    raise RuntimeError("You need to write this part!")

def spline(T):
    '''
    timeaxis, h = spline(T)
    Return a continuous spline interpolator with continuous first deriviative.

    T (scalar) - the upsampling factor
    timeaxis (array, length 4*T-1) - sample indices, from -(2*T-1) through (2*T-1)
    h (array, length 4*T-1) - the cubic spline interpolation kernel
    '''
    raise RuntimeError("You need to write this part!")

def sinc(T, D):
    '''
    timeaxis, h = sinc(T, D)
    Return D samples from the center of h(t)=sin(pi*t/T) / (pi*t/T).

    T (scalar) - the upsampling factor
    D (scalar) - the duration of the returned kernel; always an odd number.
    timeaxis (array, length D) - sample indices, from -(D-1)/2 through (D-1)/2
    h (array, length 4*T-1) - the sinc interpolation kernel
    '''
    raise RuntimeError("You need to write this part!")


