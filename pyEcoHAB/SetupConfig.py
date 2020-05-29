#!/usr/bin/env python
# encoding: utf-8
from __future__ import division, absolute_import, print_function

import os
import glob
import sys
import time

if sys.version_info < (3, 0):
    from ConfigParser import RawConfigParser, NoSectionError
else:
    from configparser import RawConfigParser, NoSectionError


class SetupConfig(RawConfigParser, matplotlib.ticker.Formatter):
    def __init__(self, path, fname=None):    
        RawConfigParser.__init__(self)
        self.path = path               
        if fname is None:
            if os.path.isfile(os.path.join(self.path, 'setup.txt')):
                self.fname = 'setup.txt'
            else:
                fnames = glob.glob(os.path.join(self.path, "setup*txt"))
                if len(fnames):
                    self.fname = fnames[0]
                else:
                    self.fname = None
        else:                  
            self.fname = fname
        if self.fname is not None:
            self.read(os.path.join(self.path, self.fname))
        else:
            pass


