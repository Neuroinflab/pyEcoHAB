#!/usr/bin/env python
# encoding: utf-8
from __future__ import division, absolute_import, print_function

import os
import glob
import sys

from collections import OrderedDict
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
                    
                    fnames = glob.glob(os.path.join(path, "setup*txt"))
                    if len(fnames):
                        self.fname = os.path.basename(fnames[0])
                        self.path = path
                    else:
                       raise Exception("No config found")

        full_path = os.path.join(self.path, self.fname)
        self.read(full_path)

    def get_cages(self):
        return filter(lambda x: x.startswith("cage"),
                      self.sections())

    def get_tunnels(self):
        return filter(lambda x: x.startswith("tunnel"),
                      self.sections())

    def get_cages_dict(self):
        secs = self.sections()
        cage_dict = OrderedDict()
        for sec in secs:
            if not sec.startswith("cage"):
                continue
            cage_dict[sec] = []
            for antenna_type, val in self.items(sec):
                if antenna_type.startswith("external"):
                    cage_dict[sec].append(val)
                if antenna_type.startswith("internal"):
                    continue
                else:
                    print("Unknown antenna type %s" % antenna_type)
            if not len(cage_dict[sec]):
                print("Did not register any antennas associated with %s", sec)
        if not len(list(cage_dict.keys())):
            print("Did not registered any cages in this setup")
        return cage_dict

    def get_pipes_dict(self):
        secs = self.sections()
        tunnel_dict = OrderedDict()
        for sec in secs:
            if not sec.startswith("tunnel"):
                continue
            tunnel_dict[sec] = []
            for antenna_type, val in self.items(sec):
                if antenna_type.startswith("entrance"):
                    tunnel_dict[sec].append(val)
                else:
                    print("Unknown antenna type %s" % antenna_type)
            if not len(tunnel_dict[sec]):
                print("Did not register any antennas associated with %s", sec)
        if not len(list(tunnel_dict.keys())):
            print("Did not registered any tunnels in this setup")
        return tunnel_dict

    def get_compartments_with_additional_antennas(self):
        out = []
        for sec in self.sections():
            all_items = self.items(sec)
            out += [item[1] for item in all_items if item[1].startswith("int")]
        return out
