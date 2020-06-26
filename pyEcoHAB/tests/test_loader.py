from __future__ import print_function, division, absolute_import
import os
import unittest
from datetime import date
import numpy as np
import pyEcoHAB.utils.for_loading as uf
import pyEcoHAB.utility_functions as utils
from pyEcoHAB import data_path, sample_data
from pyEcoHAB.SetupConfig import SetupConfig
from pyEcoHAB import Loader, Merger, Timeline
from pyEcoHAB import get_incohort_sociability
from pyEcoHAB import get_solitude
from pyEcoHAB import get_activity
from pyEcoHAB import get_dynamic_interactions

class TestLoader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path1 = os.path.join(data_path, "modular_1",
                                 "data_setup_additional")
        cls.dataset1 = Loader(cls.path1, visit_threshold=1.5, prefix="gugu")
        cls.setup1 = SetupConfig(cls.path1)
        cls.dataset1_standard = Loader(cls.path1)

    def test_path(self):
        self.assertEqual(self.path1, self.dataset1.path)

    def test_visit_threshold(self):
        self.assertEqual(self.dataset1.visit_threshold, 1.5)
    
    def test_prefix(self):
        self.assertEqual(self.dataset1_standard.prefix, "WT_F_test_1_")

    def test_prefix_2(self):
        self.assertEqual(self.dataset1.prefix, "gugu")


class TestMerger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "BALB_VPA_data_cohort_1_divided")
        cls.path1 = os.path.join(path, "setup_1")
        cls.path2 = os.path.join(path, "setup_2")
        cls.data1 = Loader(cls.path1)
        cls.data2 = Loader(cls.path2)
        cls.res_dir = os.path.join(path, "results")
        config = os.path.join(path, "experiment_setup.txt")
        cls.data = Merger(config, cls.res_dir, ecohab_1=cls.data1,
                          ecohab_2=cls.data2)
        cls.original_data = Loader(sample_data)

    def test_1(self):
        self.assertEqual(self.data.res_dir,
                         "%s_%s"%(self.res_dir,
                                  date.today().strftime("%d.%m.%y")))

    def test_cages(self):
        out = sorted(["cage A", "cage B", "cage C", "cage D"])
        self.assertEqual(sorted(self.data.cages), out)

    def test_incohort_sociability_1(self):
        config = Timeline(sample_data)
        out_1 = get_incohort_sociability(self.data, config, 3600)
        out_2 = get_incohort_sociability(self.original_data, config, 3600)
        self.assertEqual(out_1, out_2)

    def test_incohort_sociability_2(self):
        config = Timeline(sample_data)
        out_1 = get_incohort_sociability(self.data, config, 24*3600)
        out_2 = get_incohort_sociability(self.original_data, config, 24*3600)
        self.assertEqual(out_1, out_2)

    def test_solitude(self):
        config = Timeline(sample_data)
        out_1 = get_solitude(self.data, config)
        out_2 = get_solitude(self.original_data, config)
        self.assertEqual(out_1, out_2)

    def test_activity_1(self):
        config = Timeline(sample_data)
        out_1 =get_activity(self.data, config, 3600)
        out_2 =get_activity(self.original_data, config, 3600)
        self.assertEqual(out_1, out_2)

    def test_activity_2(self):
        config = Timeline(sample_data)
        out_1 = get_activity(self.data, config, 24*3600)
        out_2 = get_activity(self.original_data, config, 24*3600)
        self.assertEqual(out_1, out_2)

    def get_dynamic_interactions(self):
        config = Timeline(sample_data)
        out_1 = get_dynamic_interactions(self.data, config, N=1, seed=1)
        out_2 = get_dynamic_interactions(self.original_data, config, N=1,
                                         seed=1)
        self.assertEqual(out_1, out_2)
            
if __name__ == '__main__':
    unittest.main()

