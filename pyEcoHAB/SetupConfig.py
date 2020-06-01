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


class SetupConfig(RawConfigParser):
    STANDARD_CONFIG = {
        "cage A": {"external_antenna1": 1,
                   "external_antenna2": 8},
        "cage B": {"external_antenna1": 2,
                   "external_antenna2": 3},
        "cage C": {"external_antenna1": 4,
                   "external_antenna2": 5},
        "cage D": {"external_antenna1": 6,
                   "external_antenna2": 7},
        "pipe 1": {"antenna1": 1,
                   "antenna2": 2},
        "pipe 2": {"antenna1": 3,
                   "antenna2": 4},
        "pipe 3": {"antenna1": 5,
                   "antenna2": 6},
        "pipe 4": {"antenna1": 7,
                   "antenna2": 8},

    }
    def __init__(self, path=None, fname=None):    
        RawConfigParser.__init__(self)
        if path is None:
            self.read_dict(self.STANDARD_CONFIG)
            return
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
            self.read_dict(self.STANDARD_CONFIG)



