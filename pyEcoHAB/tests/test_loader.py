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
        cls.path2 = os.path.join(data_path, "weird_very_short")
        cls.dataset2 = Loader(cls.path2, visit_threshold=2)
        cls.setup3 = SetupConfig(path=data_path, fname="setup_short.txt")
        cls.dataset3 = Loader(cls.path1, visit_threshold=1.5,
                              setup_config=cls.setup3, remove_antennas=["8"])
        cls.path_empty = os.path.join(data_path, "empty")

    def test_load_empty(self):
        self.assertRaises(Exception, Loader, self.path_empty)

    def test_path(self):
        self.assertEqual(self.path1, self.dataset1.path)

    def test_visit_threshold(self):
        self.assertEqual(self.dataset1.visit_threshold, 1.5)

    def test_prefix(self):
        self.assertEqual(self.dataset1_standard.prefix, "WT_F_test_1_")

    def test_prefix_2(self):
        self.assertEqual(self.dataset1.prefix, "gugu")

    def test_visits_1(self):
        out = self.dataset2.get_visits()
        out2 = self.dataset2.visits.data
        self.assertEqual(len(out), len(out2))

    def test_visits_cage_A(self):
        out = self.dataset2.get_visits(cage="cage A")
        self.assertEqual(len(out), 1)

    def test_visits_cage_AA(self):
        out = self.dataset2.get_visits(cage="cage A")
        self.assertEqual(out[0].address, "cage A")

    def test_visits_cage_A_address(self):
        out = self.dataset2.get_visits(cage="cage A")
        self.assertEqual(out[0].address, "cage A")

    def test_visits_cage_A_duration(self):
        out = self.dataset2.get_visits(cage="cage A")
        self.assertTrue(np.isclose(out[0].duration,
                                   (out[0].t_end - out[0].t_start)))

    def test_visits_cage_A2(self):
        out = self.dataset2.get_visits(cage="cage A", t_end=1286708960.687)
        self.assertEqual(out, [])

    def test_visits_cage_A3(self):
        out = self.dataset2.get_visits(cage="cage A", t_start=1286708960.687)
        self.assertEqual(len(out), 1)

    def test_visits_cage_no_mice(self):
        out = self.dataset2.get_visits(mice="mouse 2")
        self.assertEqual([], out)

    def test_visits_different_cages(self):
        out = self.dataset2.get_visits(t_start=1286708669.65,
                                       t_end=1286708768.349)
        self.assertEqual(len(out), 4)

    def test_visits_different_cages_2(self):
        out = self.dataset2.get_visits(t_start=1286708669.65,
                                       t_end=1286708768.349, cage="cage C")
        self.assertEqual(len(out), 2)

    def test_visits_different_cages_3(self):
        out = self.dataset2.get_visits(t_start=1286708669.65,
                                       t_end=1286708768.349, cage="cage D")
        self.assertEqual(len(out), 2)

    def test_visits_internal_antennas_mouse1(self):
        out = self.dataset1.get_visits("mouse_1")
        out2 = self.dataset3.get_visits("mouse_1")
        self.assertEqual(out, out2)

    def test_visits_internal_antennas_mouse2(self):
        out = self.dataset1.get_visits("mouse_2")
        out2 = self.dataset3.get_visits("mouse_2")
        self.assertEqual(len(out)-1, len(out2))


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
        cls.data = Merger(config, cls.res_dir, cls.data1, cls.data2)
        cls.original_data = Loader(sample_data)

    def test_1(self):
        self.assertEqual(self.data.res_dir,
                         "%s_%s" % (self.res_dir,
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
        out_1 = get_activity(self.data, config, 3600)
        out_2 = get_activity(self.original_data, config, 3600)
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
