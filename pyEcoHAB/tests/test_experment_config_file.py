import unittest
import os

from pyEcoHAB import data_path
from pyEcoHAB import ExperimentConfigFile

class TestnoDST(unittest.TestCase):
    def test_1(self):
        path = os.path.join(data_path, "time_change")
        config = ExperimentConfigFile(path)
        times = config.gettime("5 light")
        self.assertEqual(times[1]-times[0], 12*3600)


if __name__ == '__main__':
    unittest.main()
