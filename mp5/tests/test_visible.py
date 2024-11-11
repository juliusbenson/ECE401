import unittest, h5py, submitted
from gradescope_utils.autograder_utils.decorators import weight
import numpy as np

# TestSequence
class TestStep(unittest.TestCase):
    def setUp(self):
        self.h5 = h5py.File('solutions.hdf5','r')

    @weight(7.5)
    def test_mri(self):
        with h5py.File('data.hdf5','r')  as f:
            mri_dft = f['mri_dft'][:]
        p = self.h5['downsample_and_shift_params'][:]
        mri = submitted.downsample_and_shift_dft2(mri_dft, p[0], p[1], p[2])
        e = np.sum(np.abs(mri-self.h5['mri']))/np.sum(np.abs(self.h5['mri']))
        self.assertTrue(e < 0.04, 'downsample_and_shift_dft2 wrong by more than 4% (visible case)')

    @weight(7.5)
    def test_cleaned_image(self):
        with h5py.File('data.hdf5','r')  as f:
            noisy_image = f['noisy_image'][:]
        p = self.h5['dft_filter_params'][:]
        cleaned_image = submitted.dft_filter(noisy_image, p[0], p[1], p[2], p[3])
        e = np.sum(np.abs(cleaned_image-self.h5['cleaned_image']))/np.sum(np.abs(self.h5['cleaned_image']))
        self.assertTrue(e < 0.04, 'dft_filter wrong by more than 4% (visible case)')

    @weight(7.5)
    def test_transitioned_image(self):
        with h5py.File('data.hdf5','r')  as f:
            noisy_image = f['noisy_image'][:]
        p = self.h5['dft_filter_params'][:]
        transitioned_image=submitted.transitioned_filter(noisy_image,p[0],p[1],p[2],p[3])
        e = np.sum(np.abs(transitioned_image-self.h5['transitioned_image']))/np.sum(np.abs(self.h5['transitioned_image']))
        self.assertTrue(e < 0.04, 'transitioned_filter wrong by more than 4% (visible case)')

    @weight(7.5)
    def test_zero_pad(self):
        h = self.h5['h_orig'][:]
        x = self.h5['x_orig'][:]
        hp, xp =submitted.zero_pad(h, x)
        e = np.sum(np.abs(hp-self.h5['hp']))/np.sum(np.abs(self.h5['hp']))
        self.assertTrue(e < 0.04, 'hp wrong by more than 4% (visible case)')
        e = np.sum(np.abs(xp-self.h5['xp']))/np.sum(np.abs(self.h5['xp']))
        self.assertTrue(e < 0.04, 'xp wrong by more than 4% (visible case)')

    @weight(7.5)
    def test_overlap_add(self):
        h = self.h5['h_orig'][:]
        long_signal = self.h5['long_signal'][:]
        M = 2048 + 1 - len(h)
        import time
        start_time_convolve = time.time()
        y_convolve = np.convolve(long_signal, h)
        end_time_convolve = time.time()
        complexity_convolve = end_time_convolve - start_time_convolve
        start_time_ola = time.time()
        y_ola = submitted.overlap_add(h, long_signal, M)
        end_time_ola = time.time()
        complexity_ola = end_time_ola - start_time_ola

        e = np.sum(np.abs(y_convolve-y_ola[:len(y_convolve)]))/np.sum(np.abs(y_convolve))
        self.assertTrue(e < 0.04, 'overlap_add output is wrong by more than 4% (visible case)')

        self.assertLess(complexity_ola, complexity_convolve, "overlap_add took %gsec, which doesn't beat convolve's %gsec"%(complexity_ola,complexity_convolve))
        
