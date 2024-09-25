'''
This is the module you'll submit to the autograder.

There are several function definitions, here, that raise RuntimeErrors.  You should replace
each "raise RuntimeError" line with a line that performs the function specified in the
function's docstring.
'''

import numpy as np
import cmath

def sinusoid(phasor:complex, frequency:float, duration:float, samplerate:float) -> tuple[np.array, np.array]:
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
    return np.array(timeaxis), np.array(signal)
    
def compute_aliasing(frequencies:list[float], phasors:list[complex], samplerates:list[float]) -> tuple[np.array, np.array]:
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

def fourier_analysis(signal : list[float], number_of_coefficients : int) -> np.array:
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

def interpolate(lowrate_signal : list[float], T : float, kernel_timeaxis : list[float], kernel : list[float]) -> np.array:
# def interpolate(lowrate_signal, T, kernel_timeaxis, kernel):
    '''
    highrate_signal = interpolate(lowrate_signal, T, kernel_timeaxis, kernel)
    Use lowrate-to-highrate conversion to simulate discrete-to-continuous conversion.

    lowrate_signal (length-N array) - the lowrate signal
    T (scalar) - ratio of highrate/lowrate, i.e., number of output samples per input sample
    kernel_timeaxis (array) - sample times of the kernel, at the highrate
    kernel (array) - the interpolation kernel.  length(kernel)==length(kernel_timeaxis).
    highrate_signal (length-N*T array) - the highrate signal

    Note: in order to keep the output to only N*T samples, use modulo arithmetic for the 
    interpolation, e.g.,
    highrate_signal[np.mod(kernel_timeaxis+n*T, N*T)] += kernel * lowrate_signal[n]
    '''

    N = len(lowrate_signal)
    highrate_signal = np.zeros(N*T)

    for n in range(N):
        kRange = np.mod(np.add(kernel_timeaxis, n*T), N*T)      # the + operator isn't overloaded to work with arrays like np.add is
        toPut  = np.take(highrate_signal, kRange, mode='wrap')  # python stock indexing doesn't support modular indexing
        toPut += np.multiply(kernel, lowrate_signal[n])         # the * operator isn't overloaded to work with arrays like np.multiply is

        np.put(highrate_signal, kRange, toPut, mode='wrap')

    return highrate_signal

def rectangle(T : int) -> tuple[np.array, np.array]:
    '''
    timeaxis, h = rectangle(T)
    Return a rectangle function of length T.

    T (scalar) - length, in samples
    timeaxis (length-T array) - sample indices, from 0 to T-1, corresponding to h
    h (length-T array) - the rectangle function
    '''

    timeaxis = np.arange(T)
    h = np.ones(T)
    return timeaxis, h

def triangle(T : int) -> tuple[np.array, np.array]:
    '''
    timeaxis, h = triangle(T)
    Return a triangle function of length 2*T-1.

    T (scalar) - length of each side of the triangle, in samples
    timeaxis (array, length 2*T-1) - sample indices, from -(T-1) through (T-1)
    h (array, length 2*T-1) - the triangle function, 1 - abs(timeaxis)/T
    '''
    timeaxis = np.arange(-T+1,T)
    h = np.zeros(2*T-1)

    for n,t in enumerate(timeaxis):
        h[n] = 1 - np.abs(t) / T

    return timeaxis, h

def spline(T : int) -> tuple[np.array, np.array]:
    '''
    timeaxis, h = spline(T)
    Return a continuous spline interpolator with continuous first deriviative.

    T (scalar) - the upsampling factor
    timeaxis (array, length 4*T-1) - sample indices, from -(2*T-1) through (2*T-1)
    h (array, length 4*T-1) - the cubic spline interpolation kernel
    '''

    timeaxis = np.arange(-(2*T-1),(2*T))
    h = np.zeros(len(timeaxis))
    for n,t in enumerate(timeaxis):
        if np.abs(t) <= T:
            h[n] = 1 - (3/2)*(np.abs(t)/T)**2 + (1/2)*(np.abs(t)/T)**3
        elif (T <= np.abs(t)) or (np.abs(t) <= 2*T):
            h[n] = -(3/2)*(((np.abs(t)-2*T)/T)**2)*((np.abs(t)-T)/T)
        else:
            h[n] = 0
    return timeaxis, h

def sinc(T : int, D : int) -> tuple[np.array, np.array]:
    '''
    timeaxis, h = sinc(T, D)
    Return D samples from the center of h(t)=sin(pi*t/T) / (pi*t/T).

    T (scalar) - the upsampling factor
    D (scalar) - the duration of the returned kernel; always an odd number.
    timeaxis (array, length D) - sample indices, from -(D-1)/2 through (D-1)/2
    h (array, length 4*T-1) - the sinc interpolation kernel
    '''
    timeaxis = np.arange(-(D-1)//2,((D-1)//2)+1)
    h = np.zeros(len(timeaxis))
    for n,t in enumerate(timeaxis):
        if t == 0:
            h[n] = 1
        else:
            h[n] = (np.sin(np.pi*t/T))/(np.pi*t/T)
    return timeaxis, h


