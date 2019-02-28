#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
from EcoHAB import tube_dominance as tubed
import unittest


class TestCheckTwoMice(unittest.TestCase):
    def TestLaboriusPushing_mouse2_pushing_mouse1(self):
        m1_antennas = [5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 7]
        m1_readouts  = [705.074,#5
                        708.091,#6
                        710.577,#5
                        744.813,#5
                        746.72,#5
                        758.851,#5
                        802.624,#5
                        808.095,#6
                        809.252,#6
                        809.564,#6
                        813.675]#7

       
        m2_antennas = [6, 6, 6, 6, 6, 5, 5, 5, 6, 5, 4, 4]
        m2_readouts =  [751.348,#6
                        753.771,#6
                        755.115,#6
                        755.865,#6
                        764.666,#6
                        768.856,#5
                        769.184,#5
                        794.072,#5
                        796.386,#6
                        801.342,#5
                        807.86,#4
                        814.3]#4
        out = tubed.check_mouse1_pushing_out_mouse2(m1_antennas,
                                                    m1_times,
                                                    m2_antennas,
                                                    m2_times)
        self.asserFalse(out)
    
 def TestLaboriusPushing_mouse1_pushing_mouse2(self):
        m2_antennas = [5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 7]
        m2_readouts  = [705.074,#5
                        708.091,#6
                        710.577,#5
                        744.813,#5
                        746.72,#5
                        758.851,#5
                        802.624,#5
                        808.095,#6
                        809.252,#6
                        809.564,#6
                        813.675]#7

       
        m1_antennas = [6, 6, 6, 6, 6, 5, 5, 5, 6, 5, 4, 4]
        m1_readouts =  [751.348,#6
                        753.771,#6
                        755.115,#6
                        755.865,#6
                        764.666,#6
                        768.856,#5
                        769.184,#5
                        794.072,#5
                        796.386,#6
                        801.342,#5
                        807.86,#4
                        814.3]#4
        out = tubed.check_mouse1_pushing_out_mouse2(m1_antennas,
                                                    m1_times,
                                                    m2_antennas,
                                                    m2_times)
        self.asserTrue(out)
    
