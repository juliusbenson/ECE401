import os, h5py, wave, struct
import numpy as  np

###############################################################################
# TODO: here are the functions that you need to write
def downsample_and_shift_dft2(oversampled_dft, downsampling_factor, row_shift, col_shift):
    '''
    Input: 
      oversampled_dft [M1,M2] - a 2d array containing the oversampled DFT of a grayscale image
      downsampling_factor (scalar) - the factor by which the DFT image is oversampled
      row_shift (scalar)  - the number of rows that the image should be shifted
      col_shift (scalar) - the number of columns that the image should be shifted
    Output: 
      image [M1/downsampling_factor, M2/downsampling_factor] - the real part of the inverse DFT
      of the valid frequency samples, shifted by the specified numbers of rows and columns.
    '''
    downsampled = oversampled_dft[::downsampling_factor,::downsampling_factor]
    (N1,N2) = downsampled.shape
    X_corrected = np.zeros_like(downsampled)
    for k1,row in enumerate(downsampled):
        for k2,_ in enumerate(row):
            if (0 <= k1) and (k1 < (N1/2)): # top
                rowshift = np.exp(-2j*np.pi*k1*row_shift/N1)
            else: # bottom
                rowshift = np.exp(-2j*np.pi*(k1-N1)*row_shift/N1)

            if (0 <= k2) and (k2 < (N2/2)): # left
                colshift = np.exp(-2j*np.pi*k2*col_shift/N2)
            else: # right
                colshift = np.exp(-2j*np.pi*(k2-N2)*col_shift/N2)

            X_corrected[k1,k2] = downsampled[k1,k2]*rowshift*colshift

    return np.real(np.fft.ifft2(X_corrected))
    # return np.real(X_corrected)
    # raise RuntimeError("You need to write this!")
    
def dft_filter(noisy_image, min1, max1, min2, max2):
    '''
    Input: 
      noisy_image [N1,N2] - an image with narrowband noises
      min1, max1 (scalars) - zero out all rows of the DFT min1 <= k1 < max1, likewise  for N1-k1
      min2, max2 (scalars) - zero out all cols  of the DFT min2 <= k2 < max2, likewise for N2-k2
    Outut:
      cleaned_image [N1,N2] - image with the corrupted bands removed.
      Be sure to take the real part of the inverse DFT, and then truncate
      so that 0 <= cleaned_image[n1,n2,color] <= 1 for all n1,n2,color.
    '''
    (N1,N2,_) = noisy_image.shape
    noisy_dft = np.fft.fft2(noisy_image, axes=(0,1))
    X_corrected = np.zeros_like(noisy_dft)
    for k1,row in enumerate(noisy_dft):
        for k2,_ in enumerate(row):
            if (min1 <= k1) and (k1 < max1):
                X_corrected[k1,k2] = 0
            elif (min1 <= (N1 - k1)) and ((N1 - k1) < max1):
                X_corrected[k1,k2] = 0
            elif (min2 <= k2) and (k2 < max2):
                X_corrected[k1,k2] = 0
            elif (min2 <= (N2 - k2)) and ((N2 - k2) < max2):
                X_corrected[k1,k2] = 0
            else:
                X_corrected[k1,k2] = noisy_dft[k1,k2]

    return np.real(np.fft.ifft2(X_corrected,axes=(0,1)))
    # raise  RuntimeError("You need to write this!")
    
def transitioned_filter(noisy_image, min1, max1, min2, max2):
    '''
    Input: 
      noisy_image [N1,N2] - an image with narrowband noises
      min1, max1 (scalars) - zero out all rows of the DFT min1 <= k1 < max1, likewise  for N1-k1
      min2, max2 (scalars) - zero out all cols  of the DFT min2 <= k2 < max2, likewise for N2-k2
    Outut:
      cleaned_image [N1,N2] - image with the corrupted bands removed.
      Be sure to take the real part of the inverse DFT, and then truncate
      so that 0 <= cleaned_image[n1,n2,color] <= 1 for all n1,n2,color.

    Transition band:
      the bands k1=min1-1, k1=max1, k2=min2-1, and k2=max2 should be set to half of their
      original values, 0.5*X[k1,k2].
    '''
    (N1,N2,_) = noisy_image.shape
    noisy_dft = np.fft.fft2(noisy_image, axes=(0,1))
    X_corrected = np.zeros_like(noisy_dft)
    for k1,row in enumerate(noisy_dft):
        for k2,_ in enumerate(row):
            if (min1 <= k1) and (k1 < max1):
                X_corrected[k1,k2] = 0
            elif (min1 <= (N1 - k1)) and ((N1 - k1) < max1):
                X_corrected[k1,k2] = 0
            elif (min2 <= k2) and (k2 < max2):
                X_corrected[k1,k2] = 0
            elif (min2 <= (N2 - k2)) and ((N2 - k2) < max2):
                X_corrected[k1,k2] = 0
            elif (min1-1 <= k1) and (k1 < max1):
                X_corrected[k1,k2] = noisy_dft[k1,k2]/2
            elif (min2-1 <= k2) and (k2 < max2):
                X_corrected[k1,k2] = noisy_dft[k1,k2]/2
            elif (min1-1 <= (N1 - k1)) and ((N1 - k1) < max1):
                X_corrected[k1,k2] = noisy_dft[k1,k2]/2
            elif (min2-1 <= (N2 - k2)) and ((N2 - k2) < max2):
                X_corrected[k1,k2] = noisy_dft[k1,k2]/2
            else:
                X_corrected[k1,k2] = noisy_dft[k1,k2]

    return np.real(np.fft.ifft2(X_corrected,axes=(0,1)))
    # raise  RuntimeError("You need to write this!")

def zero_pad(h, x):
    '''
    (hp,xp) = zero_pad(h,x)
    Input:
      h [L] - a length-L impulse response array
      x [M] - a length-M signal array
    Return: 
      hp [N] - the same h, but zero-padded to a length of N=L+M-1
      xp [N] - the same x, but zero-padded to a length of N=L+M-1
    '''
    L = len(h)
    M = len(x)
    N = L + M - 1

    hp = np.zeros(N)
    xp = np.zeros(N)

    hp[:L]=h[:]
    xp[:M]=x[:]

    return hp, xp
    # raise  RuntimeError("You need to write this!")

def overlap_add(h, x, M):
    '''
    y = overlap_add(h, x, M)
    Input:
      h [L] - a length-L impulse response array
      x - a very long array containing the signal x
      M (scalar) - the length of frames into which x should be chopped
    Output:
      y - the result of filtering x by h using overlap-add

    First, zero-pad h to a length of N samples, where N=L+M-1.
    Second, compute its DFT, H.
    Third, chop x into frames of length M.  
      There should be NF of these; the last one might have fewer than M samples.
    Fourth, prepare the output array, y, of length equal to NF*M + L - 1.
    Fifth, for each frame of M samples in x:
     - zero-pad to a length of N
     - compute its DFT
     - multiply its DFT by H
     - inverse DFT
     - add the inverse DFT to the correct place in the y array
    Return y
    '''
    L = len(h)
    N = L+M-1
    NF = (int)(np.ceil(len(x) / M))

    # Zero-pad x
    if len(x) < NF*M:
        xp = np.concatenate((x,np.zeros(NF*M-len(x)))) # This will make things easier later

    # First, zero-pad h to a length of N samples, where N=L+M-1.
    hp = np.zeros(N)
    hp[:L]=h[:]

    # Second, compute its DFT, H.
    H = np.fft.fft(hp)

    # Third, chop x into frames of length M.
    frames = np.split(xp,NF)

    # Fourth, prepare the output array, y, of length equal to NF*M + L - 1.
    y = np.zeros(NF*M + L - 1)

    # Fifth, for each frame of M samples in x:
    for i,frame in enumerate(frames):
    #  - zero-pad to a length of N
        fp = np.zeros(N)
        fp[:len(frame)]=frame[:]

    #  - compute its DFT
        F = np.fft.fft(fp)

    #  - multiply its DFT by H
        FH = np.multiply(F,H)

    #  - inverse DFT
        fh = np.fft.ifft(FH)
        # print(f'len(fh) {len(fh)}')

    #  - add the inverse DFT to the correct place in the y array
        y[i*M:i*M+N]=y[i*M:i*M+N]+fh[:]
    return y

    # raise  RuntimeError("You need to write this!")
