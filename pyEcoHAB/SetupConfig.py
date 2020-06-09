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
                    
                    fnames = glob.glob(os.path.join(path, "setup*txt"))
                    if len(fnames):
                        self.fname = os.path.basename(fnames[0])
                        self.path = path
                    else:
                       raise Exception("No config found")

        full_path = os.path.join(self.path, self.fname)
        self.read(full_path)

        self.cages = self.get_cages()
        self.tunnels = self.get_tunnels()
        self.cages_dict = self.get_cages_dict()
        self.tunnels_dict = self.get_tunnels_dict()
        self.same_tunnel = self.get_same_tunnel()
        self.same_address = self.get_same_address()
        self.opposite_tunnel = self.get_opposite_tunnel_dict()
        self.address = self.get_cage_address_dict()
        self.address_non_adjacent = self.get_address_non_adjacent_dict()
        self.address_surrounding = self.get_surrounding_dict()
        self.directions = self.get_directions_dict()

    def get_cages(self):
        return sorted(filter(lambda x: x.startswith("cage"),
                      self.sections()))

    def get_tunnels(self):
        return sorted(filter(lambda x: x.startswith("tunnel"),
                      self.sections()))

    def get_cages_dict(self):
        cage_dict = {}
        cages = self.get_cages()
        for sec in cages:
            cage_dict[sec] = []
            for antenna_type, val in self.items(sec):
                if antenna_type.startswith("entrance"):
                    cage_dict[sec].append(val)
                elif antenna_type.startswith("internal"):
                    continue
                else:
                    print("Unknown antenna type %s" % antenna_type)
            if not len(cage_dict[sec]):
                print("Did not register any antennas associated with %s", sec)
        if not len(list(cage_dict.keys())):
            print("Did not registered any cages in this setup")
        return cage_dict

    def get_tunnels_dict(self):
        tunnel_dict = {}
        tunnels = self.get_tunnels()
        for sec in tunnels:
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

    @property
    def internal_antennas(self):
        out = []
        for sec in self.sections():
            all_items = self.items(sec)
            out += [sec for item in all_items if item[0].startswith("int")]
        return out

    def get_same_tunnel(self):
        tunnel_dict = self.get_tunnels_dict()
        out = {}
        for tunnel, value in tunnel_dict.items():
            for antenna in value:
                out[antenna] = value
        return out

    def get_same_address(self):
        cage_dict = self.get_cages_dict()
        out = {}
        for cage, value in cage_dict.items():
            for antenna in value:
                out[antenna] = value
        return out

    @property
    def entrance_antennas(self):
        out = []
        for sec in self.sections():
            for key, value in self.items(sec):
                if key.startswith("entrance") and value not in out:
                    out.append(value)
        return out

    def other_tunnel_antenna(self, antenna):
        tunnel_antennas = self.same_tunnel[antenna][:]
        idx = tunnel_antennas.index(antenna)
        tunnel_antennas.pop(idx)
        return tunnel_antennas

    def other_cage_antenna(self, antenna):
        cage_antennas = self.same_address[antenna][:]
        idx = cage_antennas.index(antenna)
        cage_antennas.pop(idx)

        return cage_antennas

    def next_tunnel_antennas(self, antenna):
        out = []
        same_pipe = self.same_tunnel[antenna]
        for ant in same_pipe:
            same_cage_a = self.other_cage_antenna(ant)
            for a_2 in same_cage_a:
                other_pipe = self.same_tunnel[a_2]
                if other_pipe != same_pipe:
                    out.extend(other_pipe)
        return sorted(out)

    def get_opposite_tunnel_dict(self):
        # distance equal two
        all_antennas = self.entrance_antennas
        same_cages = self.same_address
        same_pipe = self.same_tunnel
        out = {}
        for a_1 in all_antennas:
            same_cage_antennas = self.other_cage_antenna(a_1)

            for a_2 in same_cage_antennas:
                pipe_next = self.other_tunnel_antenna(a_2)

                for a_3 in pipe_next:
                    cage_plus_2 = self.other_cage_antenna(a_3)

                    for a_4 in cage_plus_2:
                        tunnel_antennas = same_pipe[a_4]
                        next_tunnel_antennas = self.next_tunnel_antennas(a_4)
                        if a_1 not in tunnel_antennas and a_1 not in next_tunnel_antennas:
                            if a_1 not in out:
                                out[a_1] = []
                            for ant in tunnel_antennas:
                                if ant not in out[a_1]:
                                    out[a_1].extend(ant)

        return out

    def get_cage_address_dict(self):
        out = {}
        for sec in self.cages:
            for antenna_type, antenna in self.items(sec):
                if antenna_type.startswith("entrance"):
                    if antenna in out:
                        raise Exception("%s was specified as %s twice"%(antenna_type,
                                                                        antenna))
                    else:
                        out[antenna] = sec
        return out

    def get_address_non_adjacent_dict(self):
        all_antennas = self.entrance_antennas
        out = {}
        cage_dict = self.get_cage_address_dict()
        for antenna in all_antennas:
            pipe_next = self.other_tunnel_antenna(antenna)
            try:
                cage_adjacent = cage_dict[pipe_next[0]]
                out[antenna] = cage_adjacent
            except:
                pass
        return out

    def get_surrounding_dict(self):
        all_antennas = self.entrance_antennas
        out = {}
        cage_dict = self.get_cage_address_dict()

        for antenna in all_antennas:
            pipe_next = self.other_tunnel_antenna(antenna)
            try:
                cage_adjacent_antennas = self.other_cage_antenna(pipe_next[0])
            except:
                continue
            for caa in cage_adjacent_antennas:
                key = (min(antenna, caa), max(antenna, caa))
                if key not in out:
                    out[key] = cage_dict[caa]
        return out

    def get_directions_dict(self):
        out = []
        for tunnel in self.tunnels:
            vals = [item[1] for item in self.items(tunnel) if item[0].startswith("entra")]
            if len(vals) > 2:
                raise Exception("There are more than 2 antennas at the entrances to %s" % tunnel)
            out += [vals[0]+vals[1], vals[1]+vals[0]]
        return sorted(out)

