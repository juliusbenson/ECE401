import numpy as np

def todo_spectrum(x):
    '''
    Input:
    x (N) - a waveform
    Output:
    spec (N) - magnitude spectrum at N frequencies between 0 and fs, not including fs.
    '''
    k = np.fft.fft(x)
    return np.abs(k)
    raise RuntimeError('You need to write this!')

def todo_findpeak(spec, fs, flo, fhi):
    '''
    Input:
    spec (N) - magnitude spectrum at N frequencies between 0 and fs, not including fs.
    fs (scalar)  - sampling frequency
    flo (scalar) -  low end of the frequency range to search
    fhi (scalar) - high end of the frequency range to search
    Output:
    fpeak (scalar) - frequency of the highest-magnitude spectral sample, in Hertz
    '''
    f = np.linspace(0,fs,len(spec),endpoint=False)
    subspec = spec[(flo<f) & (f<fhi)]
    arg = np.argmax(subspec)
    fstart = len([a for a in f if a <= flo])
    return f[arg + fstart]
    raise RuntimeError('You need to write this!')

def todo_zeros(freqs, fs):
    '''
    Input:
    freqs (nfreqs) - an array of frequencies you want zeroed out, in Hertz
    fs (scalar) - sampling frequency, in Hertz
    Output:
    r (2,nfreqs) - an array of complex zeros on the unit circle, in complex conjugate pairs
    '''
    r_0k = np.exp(np.multiply(freqs,2j*np.pi/fs))
    return r_0k,np.conjugate(r_0k)
    raise RuntimeError('You need to write this!')

def todo_poles(r,BW,fs):
    '''
    Input: 
    r (2,nfreqs) - an array of complex zeros on the unit circle, in complex conjugate pairs
    BW (scalar) - desired bandwidth, in Hertz
    fs (scalar) - sampling frequency, in Hertz
    Output:
    p (2,nfreqs) - an array of complex poles with bandwidth BW, in complex conjugate pairs
    '''
    a = np.exp(-BW*np.pi/fs)
    return np.multiply(a,r)
    raise RuntimeError('You need to write this!')

def todo_coefficients(r):
    '''
    Input: 
    r (2,nfreqs) - an array of complex roots, in complex conjugate pairs
    Output:
    b (3,nfreqs) - an array of second-order polynomial coefficients, one per complex root pair
    '''
    (roots,conjs) = r
    nfreqs = len(roots)
    b = np.zeros((nfreqs,3))
    for i,root in enumerate(roots):
        conj = conjs[i]
        b[i] = [1,-root-conj,(-root)*(-conj)]
    return b.T
    raise RuntimeError('You need to write this!')

def todo_freqresponse(a, b, N, version='right'):
    '''
    Input: 
    a (3) - feedback coefficients.  You may assume a[0]=1.
    b (3) - feedforward coefficients.  You may assume b[0]=1.
    N (scalar) - number of samples of the frequency response to compute
    Output: 
    omega (N) - frequencies linearly spaced between 0 and 2pi, not including 2pi.
    H (N) - B(e^{jw})/A(e^{jw}) evaluated at the frequencies in the vector omega.
    '''
    omega = np.linspace(0,2*np.pi,N,endpoint=False)
    H = np.ones(N,dtype='complex')
    for n,w in enumerate(omega):
        A = np.sum([a_ik*np.exp(1j*w*-i) for i,a_ik in enumerate(a)])
        B = np.sum([b_ik*np.exp(1j*w*-i) for i,b_ik in enumerate(b)])
        H[n] = np.divide(B,A)
    return omega, H
    raise RuntimeError('You need to write this!')
    
def todo_filter(x, a, b):
    '''
    Input: 
    a (3) - feedback coefficients.  You may assume a[0]=1.
    b (3) - feedforward coefficients.  You may assume b[0]=1.
    x (N) - input waveform
    Output: 
    y (N) - output after being filtered using B(z)/A(z)
      Assume that x[n]==0 for n<0.
      Do not generate samples of y[n] for n >= N.
    '''
    N = len(x)
    y = np.zeros(N)
    for n in range(N):
        if n > 2:
            y[n] = x[n] + b[1]*x[n-1] + b[2]*x[n-2] - a[1]*y[n-1] - a[2]*y[n-2]
        elif n > 1:
            y[n] = x[n] + b[1]*x[n-1]               - a[1]*y[n-1]
        else:
            y[n] = x[n]
    return y
    raise RuntimeError('You need to write this!')

def todo_bell_pole(sample_rate, frequency, decay_time):
    '''
    Compute the pole-pair for a damped resonator with given frequency and decay time.
    Note that bandwidth of a pole is 1/(pi * decay_time).

    @param:
    sample_rate (int) - sampling rate, in samples/second
    frequency (int) - frequency of the resonator, in Hertz
    decay_time (int) - decay time of the resonator, in seconds
    @return:
    p (complex 2-array) - poles in the Z plane
    '''
    BW = 1/np.multiply(np.pi,decay_time)
    r = todo_zeros(frequency,sample_rate)
    p = todo_poles(r,BW,sample_rate)
    return p
    raise RuntimeError('You need to write this!')

def todo_bell_simple(sample_rate, frequency, decay_time, duration):
    '''
    Generate a decaying resonant impulse response at the given frequency and decay time.
    The numerator coefficients are always b = [1,0,0]
    The denominator coefficients are 
    a = todo_coefficinets(todo_bell_pole(sample_rate,frequency,decay_time))

    Input:
    sample_rate (integer) - sampling rate, in samples/second
    frequency (nfreqs) - frequency of the resonator, in Hertz
    decay_time (nfreqs) - decay time of the resonator, ni seconds
    duration (float) - duration of x[n], in seconds
    
    Output:
    x (array of length sample_rate*duration) - impulse response of the resonator
    '''
    (root,conj) = todo_bell_pole(sample_rate,frequency,decay_time)
    a = todo_coefficients(([root],[conj]))
    a = np.concat(a)
    b = [1, 0, 0]
    samples = int(duration*sample_rate)
    d = np.zeros(samples)
    for i in range(samples):
        d[i] = np.cos(660/sample_rate)
    # d[0] = 1 # kronecker delta = discrete impulse
    x = todo_filter(d, a, b)
    return x
    raise RuntimeError('You need to write this!')

def todo_bell_multitone(sample_rate, frequencies, decay_times, duration):
    '''
    Generate the sound of a bell ringing at several frequencies simultaneously.

    Input:
    sample_rate (integer) - sampling rate, in samples/second
    frequencies (nfreqs) - frequencies of the bell, in Hertz
    decay_times (nfreqs) - times after onset (in seconds) at which each tone decays to exp(-1)
    duration (float) - duration of x[n], in seconds
    
    Output:
    x (array of length sample_rate*duration) - impulse response of the bell.
      This is just the sum of the todo_bell_simple impulse responses.
    '''
    p = todo_bell_pole(sample_rate,frequencies,decay_times)
    a = todo_coefficients(p)
    b = [1,0,0]
    samples = int(duration*sample_rate)
    d = np.zeros(samples)

    for i in range(samples):
        d[i] = np.cos(660/sample_rate)

    x = d
    for k in range(len(frequencies)):
        y = todo_filter(x, a[:,k], b)
        x = y

    return y
    raise RuntimeError('You need to write this!')
