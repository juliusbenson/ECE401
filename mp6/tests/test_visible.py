import unittest, h5py, submitted
from gradescope_utils.autograder_utils.decorators import weight
import numpy as np

# TestSequence
class TestStep(unittest.TestCase):
    def setUp(self):
        self.h5 = h5py.File('solutions.hdf5','r')
        self.fs = self.h5['fs'][0]
        self.duration = self.h5['duration'][0]

    @weight(3.75)
    def test_spectrum(self):
        spec=submitted.todo_spectrum(self.h5['x'])
        e = np.sum(np.abs(spec-self.h5['spec']))/np.sum(np.abs(self.h5['spec']))
        self.assertTrue(e < 0.04, 'todo_spec wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_findpeak(self):
        noise1 = submitted.todo_findpeak(self.h5['spec'], self.fs, 0, 100)
        noise2 = submitted.todo_findpeak(self.h5['spec'], self.fs, noise1+1, 100)
        noise3 = submitted.todo_findpeak(self.h5['spec'], self.fs, 100, 110)
        noise4 = submitted.todo_findpeak(self.h5['spec'], self.fs, 100, 150)
        freqs = np.array([noise1,noise2,noise3,noise4])    
        e = np.sum(np.abs(freqs-self.h5['freqs']))/np.sum(np.abs(self.h5['freqs']))
        self.assertLess(e, 0.04, 'todo_findpeak wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_zeros(self):
        z= submitted.todo_zeros(self.h5['freqs'], self.fs)
        e = np.sum(np.abs(z-self.h5['r']))/np.sum(np.abs(self.h5['r']))
        self.assertLess(e, 0.04, 'todo_zeros wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_poles(self):
        p= submitted.todo_poles(self.h5['r'],20,self.fs)
        e = np.sum(np.abs(p-self.h5['p']))/np.sum(np.abs(self.h5['p']))
        self.assertLess(e, 0.04, 'todo_poles wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_coefficients(self):
        a = submitted.todo_coefficients(self.h5['p'])
        b = submitted.todo_coefficients(self.h5['r'])
        e = np.sum(np.abs(a-self.h5['a']))/np.sum(np.abs(self.h5['a']))
        self.assertLess(e, 0.04, 'todo_coefficients: "a" wrong by more than 4% (visible case)')
        e = np.sum(np.abs(b-self.h5['b']))/np.sum(np.abs(self.h5['b']))
        self.assertLess(e, 0.04, 'todo_coefficients: "b" wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_freqresponse(self):
        nsamps = len(self.h5['x'])
        H = np.ones(nsamps,dtype='complex')
        for k in range(len(self.h5['freqs'])):
            omega, Hnew = submitted.todo_freqresponse(self.h5['a'][:,k], self.h5['b'][:,k], nsamps)
            H *= Hnew
        e = np.sum(np.abs(H-self.h5['H']))/np.sum(np.abs(self.h5['H']))
        e2 = np.sum(np.abs(H-self.h5['H_wrong']))/np.sum(np.abs(self.h5['H_wrong']))
        self.assertTrue((e<0.04) or (e2<0.04), 'todo_freqresponse wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_filter(self):
        x_tmp = self.h5['x'][:]
        for k in range(len(self.h5['freqs'])):
            y_tmp = submitted.todo_filter(x_tmp, self.h5['a'][:,k], self.h5['b'][:,k])
            x_tmp = y_tmp
        y = y_tmp
        e = np.sum(np.abs(y-self.h5['y']))/np.sum(np.abs(self.h5['y']))
        self.assertLess(e, 0.04, 'todo_filter wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_bell_pole(self):
        p = submitted.todo_bell_pole(self.fs, self.h5['bell_frequencies'][0],
                                      self.h5['bell_decays'][0])
        e = np.sum(np.abs(p-self.h5['bell_poles']))/np.sum(np.abs(self.h5['bell_poles']))
        self.assertLess(e, 0.04, 'todo_bell_poles wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_bell_simple(self):
        h = submitted.todo_bell_simple(self.fs, self.h5['bell_frequencies'][0],
                                       self.h5['bell_decays'][0], self.duration)
        self.assertEqual(len(self.h5['bell_simple'][:]), len(h),
                         'todo_bell_simple should have length %d'%(len(self.h5['bell_simple'])))
        e = np.sum(np.abs(h-self.h5['bell_simple']))/np.sum(np.abs(self.h5['bell_simple']))
        self.assertLess(e, 0.04, 'todo_bell_simple wrong by more than 4% (visible case)')

    @weight(3.75)
    def test_bell_multitone(self):
        h = submitted.todo_bell_multitone(self.fs, self.h5['bell_frequencies'][:],
                                          self.h5['bell_decays'][:], self.duration)
        self.assertEqual(len(self.h5['bell_multitone'][:]), len(h),
                         'todo_bell_multitone should have length %d'%(len(self.h5['bell_multitone'])))
        e = np.sum(np.abs(h-self.h5['bell_multitone']))/np.sum(np.abs(self.h5['bell_multitone']))
        self.assertLess(e, 0.04, 'todo_bell_simple wrong by more than 4% (visible case)')

