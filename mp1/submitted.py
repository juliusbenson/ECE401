'''
This is the module you'll submit to the autograder.

There are several function definitions, here, that raise RuntimeErrors.  You should replace
each "raise RuntimeError" line with a line that performs the function specified in the
function's docstring.
'''

import numpy as np
import note2f0
import time

def level_to_amplitude(levels):
    '''
    amplitudes = level_to_amplitudes(levels)

    levels - any array or list of levels, in dB
    amplitudes - a numpy array or list, of the same length, containing amplitudes
    '''

    amps = []
    for level in levels:
        amps.append(10 ** (level / 20))
    return amps

def create_phasors(amplitudes, phases):
    '''
    phasors = create_phasors(amplitudes, phases)
    
    phases - an array or list of the same length, containing phase angles, in radians
    amplitudes - an array or list of amplitudes
    phasors - the resulting phasors = amplitude*exp(j*phases)
    '''
    # raise RuntimeError("You need to write this part!")

    if len(amplitudes) == len(phases):
        phasors = []
        for i,_ in enumerate(amplitudes):
            # amplitude (radius) and phase (theta) convert to cartesian coordinates
            phase = phases[i]
            amplitude = amplitudes[i]
            re = amplitude*np.cos(phase)
            im = amplitude*np.sin(phase)
            phasor = complex(real=re,imag=im)
            phasors.append(phasor)
        return phasors
    else:
        raise RuntimeError("Expected same number of amplitudes and phases. Got%d amps and %d phases" % (len(amplitudes),len(phases),))

def synthesize_note(z, F0, Fs, d):
    '''
    x = synthesize_note(z, F0, Fs, d)
    
    z (array of length N) - an array or list of phasors, giving the amplitude and phase of each harmonic
    F0 (scalar) - the fundamental frequency, in Hertz
    Fs (scalar) - the sampling frequency, in samples/second
    d (scalar) - the duration, in seconds
    x (array of length Fs*d) - the synthesized signal
    
    This function creates one harmonic for each element of z, then adds them together to generate x.
    '''

    x = []
    samples = (int)(Fs * d)
    rs = np.absolute(z)
    ps = np.angle(z)
    for sample in range(0,samples):
        a = []
        t = sample / Fs # time at this sample
        for k,_ in enumerate(z):
            # r = np.absolute(phasor)
            # p = np.angle(phasor)
            a.append((rs[k] * np.cos(ps[k] + (2 * np.pi * (k + 1) * F0 * t))))
        x.append(np.sum(a))
    return x

def names_to_fundamentals(names):
    '''
    F0 = names_to_fundamentals(names)
    
    names - a list of names of notes, e.g., ['D#4','G#4','F4','G4','F4']
    fundamentals - a list of the corresponding fundamental frequencies, in Hertz, e.g., [311.13, 415.3, 349.23, 392.0, 349.23]
    '''

    F0 = []
    for name in names:
        F0.append(note2f0.F0[name])
    return F0

def synthesize_song(fundamentals, beats, bpm, Fs, phasors):
    '''
    song = synthesize_song(fundamentals, beats, bpm, Fs, phasors)
    
    fundamentals (array) - fundamental frequencies of the notes to be played
    beats (array) - duration of each note, in beats, e.g., [1,3/4,1/4,1,1]
    bpm (scalar) - number of beats per minute
    Fs (scalar) - number of samples per second
    phasors (list or array) - amplitude and phase of each harmonic
    song (numpy array)  - the returned song
    
    This function should first use beats and bpm to figure out the durations of the notes, in seconds.
    Then, for each note, it should:
        (1) call synthesize_note, to synthesize the note
        (2) call np.hanning, to create a hanning window, and multiply it by the note
        (3) call np.concatenate to concatenate the note  onto the song
    '''

    durations = []
    for i,beat in enumerate(beats):
        duration = 60 * beat / bpm # duration in seconds
        print(f'beat {i}: {duration}\t= 60 * {beat} / {bpm}')
        durations.append(duration)

    notes = []
    for n,d in enumerate(durations):
        t0 = time.time()
        fundamental = fundamentals[n]
        note = synthesize_note(phasors,fundamental,Fs,d)
        hanning = np.hanning(len(note))
        notes.append(note * hanning)
        t1 = time.time()
        te = t1-t0
        print(f'{d}\tsec note ({d*Fs} samples) with {len(phasors)} harmonics synthesized and appended in {te} sec')
    return np.concatenate(notes)
