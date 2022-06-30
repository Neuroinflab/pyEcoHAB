import unittest
import pyEcoHAB.utils.temporal as ut

class TestConvertIntsToTime(unittest.TestCase):
    def test_1(self):
        self.assertEqual(ut.convert_int_to_time(1), "01")

    def test_2(self):
        self.assertEqual(ut.convert_int_to_time(0), "00")

    def test_3(self):
        self.assertEqual(ut.convert_int_to_time(31), "31")

class TestFindLightBeginnig(unittest.TestCase):
    def test_1(self):
        dark_beg = "12:00"
        dark_len = 12
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "00:00")
        
    def test_2(self):
        dark_beg = "12:40"
        dark_len = 12.5
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "01:10")


    def test_3(self):
        dark_beg = "00:00"
        dark_len = 12.5
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "12:30")

if __name__ == '__main__':
    unittest.main()
