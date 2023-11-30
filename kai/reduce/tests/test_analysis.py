import unittest
from kai.reduce import analysis

class TestAlignCombo(unittest.TestCase):
    def test_alignComboFlyStar(self):
        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        root_dir = mod_pat + '/../data/test_epoch/'
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        combo_dir = epoch_dir + 'combo/'

        analysisObject = analysis.Analysis('17may21', filt='kp', rootDir=root_dir)



        return


if __name__ == '__main__':
    unittest.main()
