import unittest
from kai.reduce import lin_correction
from kai import instruments


class LinearityCorrectionTests(unittest.TestCase):
    def test_osiris_linearity(self):
        osiris = instruments.OSIRIS()

        # This should successfully run and no error should be thrown.
        lin_correction.lin_correction('foo.fits', instrument=osiris)

        return

if __name__ == '__main__':
    unittest.main()
