#!/usr/bin/env python
# encoding: utf-8
from __future__ import division, absolute_import, print_function

import os
import glob
import sys


from pyEcoHAB import data_path

if sys.version_info < (3, 0):
    from ConfigParser import RawConfigParser, NoSectionError
else:
    from configparser import RawConfigParser, NoSectionError


class SetupConfig(RawConfigParser):
    def __init__(self, path=None, fname=None):    
        RawConfigParser.__init__(self)

        if path is None:
            self.path = data_path
            self.fname = "standard_setup.txt"
        else:
            self.path = path
            if fname is not None:
                self.fname = fname
            else:
                if os.path.isfile(os.path.join(self.path, 'setup.txt')):
                    self.fname = 'setup.txt'
                else:
                    fnames = glob.glob(os.path.join(self.path, "setup*txt"))
                    if len(fnames):
                        self.fname = fnames[0]
                    else:
                        self.path = data_path
                        self.fname = "standard_setup.txt"

        full_path = os.path.join(self.path, self.fname)
        self.read(full_path)


   
