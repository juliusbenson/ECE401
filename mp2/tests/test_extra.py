import unittest, h5py, extra
from gradescope_utils.autograder_utils.decorators import weight
import numpy as np

# TestSequence
class TestStep(unittest.TestCase):
    @weight(3.75)
    def test_extra(self):
        with h5py.File('extra_solutions.hdf5','r') as f:
            lowrate_signal, coefficients = extra.rate_conversion(f['highrate_signal'][:],
                                                                 f['highrate'][0], f['lowrate'][0])
            self.assertAlmostEqual(
                np.average(np.abs(f['lowrate_signal']-lowrate_signal)), 0, places=3,
                msg="\n*** extra credit: lowrate_signal is wrong by more than 0.001 average"
            )
            self.assertAlmostEqual(
                np.average(np.abs(f['coefficients']-coefficients)), 0, places=3,
                msg="\n*** extra credit: coefficients wrong by more than 0.001 average"
            )
