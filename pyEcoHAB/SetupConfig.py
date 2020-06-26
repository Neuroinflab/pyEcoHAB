#!/usr/bin/env python
# encoding: utf-8
from __future__ import division, absolute_import, print_function

import os
import glob
import sys

from pyEcoHAB import data_path

if sys.version_info < (3, 0):
    from ConfigParser import RawConfigParser, NoSectionError, DuplicateSectionError
else:
    from configparser import RawConfigParser, NoSectionError, DuplicateSectionError

class SetupConfigMethods(RawConfigParser):
    
    def __init__(self):
        RawConfigParser.__init__(self)

    def make_definitions(self):
        self.all_antennas = self.get_all_antennas()
        self.cages_dict = self.get_cages_dict()
        self.tunnels_dict = self.get_tunnels_dict()
        self.same_tunnel = self.get_same_tunnel()
        self.same_address = self.get_same_address()
        self.opposite_tunnel = self.get_opposite_tunnel_dict()
        self.address = self.get_cage_address_dict()
        self.address_non_adjacent = self.get_address_non_adjacent_dict()
        self.address_surrounding = self.get_surrounding_dict()
        self.directions = self.get_directions_dict()

    def get_all_antennas(self):
        all_antennas = []
        for sec in self.sections():
            if sec.startswith("shared point"):
                continue
            for key, value in self.items(sec):
                if value not in all_antennas:
                    all_antennas.append(value)
        return sorted(all_antennas)

    @property
    def cages(self):
        return sorted(filter(lambda x: "cage" in x,
                      self.sections()))

    @property
    def tunnels(self):
        return sorted(filter(lambda x: "tunnel" in x,
                      self.sections()))

    def get_cages_dict(self):
        cage_dict = {}
        for sec in self.cages:
            cage_dict[sec] = []
            for antenna_type, val in self.items(sec):
                if antenna_type.startswith("entrance"):
                    cage_dict[sec].append(val)
                elif antenna_type.startswith("internal"):
                    cage_dict[sec].append(val)
                else:
                    print("Unknown antenna type %s" % antenna_type)
            if not len(cage_dict[sec]):
                print("Did not register any antennas associated with %s", sec)
        if not len(list(cage_dict.keys())):
            print("Did not registered any cages in this setup")
        return cage_dict

    def get_tunnels_dict(self):
        tunnel_dict = {}
        for sec in self.tunnels:
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
            out += [item[1] for item in all_items if item[0].startswith("int")]
        return out

    def get_same_tunnel(self):
        tunnel_dict = self.get_tunnels_dict()
        out = {}
        for tunnel, value in tunnel_dict.items():
            for antenna in value:
                out[antenna] = [val for val in value]
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

    def other_tunnel_antenna(self, new_antenna):
        antenna = new_antenna
        tunnel_antennas = self.same_tunnel[antenna][:]
        idx = tunnel_antennas.index(antenna)
        tunnel_antennas.pop(idx)
        return tunnel_antennas

    def other_cage_antenna(self, new_antenna):
        antenna = new_antenna
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
                if a_2 in self.internal_antennas:
                    continue
                other_pipe = self.same_tunnel[a_2]
                if other_pipe != same_pipe:
                    out.extend(other_pipe)
        return sorted(out)

    def _go_two_steps(self, antenna):
        out = []
        for a_2 in self.other_cage_antenna(antenna):
            if a_2 in self.internal_antennas:
                continue
            pipes_next = self.other_tunnel_antenna(a_2)
            for a_3 in pipes_next:
                if a_3 in self.internal_antennas:
                    continue
                cage_plus_2 = self.other_cage_antenna(a_3)
                for a_4 in cage_plus_2:
                    if a_4 in self.internal_antennas:
                        continue
                    tunnel_antennas = self.same_tunnel[a_4]
                    next_tunnel_antennas = self.next_tunnel_antennas(a_4)
                    if antenna not in tunnel_antennas and antenna not in next_tunnel_antennas:
                        out += tunnel_antennas
        return list(set(out))

    def get_opposite_tunnel_dict(self):
        # distance equal two
        all_antennas = self.entrance_antennas
        out = {}
        for a_1 in all_antennas:
            out_this_antenna = self._go_two_steps(a_1)
            other_tunnel_antennas = self.other_tunnel_antenna(a_1)
            out_other_tunnel_antenna = []
            for other_tunnel_antenna in other_tunnel_antennas:
                out_other_tunnel_antenna += self._go_two_steps(other_tunnel_antenna)
            if len(out_this_antenna + out_other_tunnel_antenna):
                out[a_1] = sorted(list(set(out_this_antenna + out_other_tunnel_antenna)))
        return out

    def get_cage_address_dict(self):
        out = {}
        for sec in self.cages:
            for antenna_type, antenna in self.items(sec):
                if antenna_type.startswith("entrance") or antenna_type.startswith("internal"):
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
        out = {}
        cage_dict = self.get_cage_address_dict()
        for antenna in self.entrance_antennas:
            pipe_next = self.other_tunnel_antenna(antenna)
            cage_adjacent_antennas = self.other_cage_antenna(pipe_next[0])
            for caa in cage_adjacent_antennas:
                if caa == antenna:
                    continue
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
            out += [vals[0]+"_"+vals[1], vals[1]+"_"+vals[0]]
        return sorted(out)

    def find_unused_antennas(self):
        out_l = []
        for sec in self.sections():
            ants = [item[1] for item in self.items(sec)]
            out_l.extend(ants)
        return sorted(set(self.ALL_ANTENNAS) - set(out_l))

    def get_mismatched_pairs(self):
        pairs = []
        for i in range(1, 9):
            for j in range(i+1, 9):
                pairs.append("%d %d" % (i, j))

        unused = sorted(list(self.find_unused_antennas()))
        legal = []
        for i, u in enumerate(unused):
            for u2 in self.ALL_ANTENNAS:
                if u == u2:
                    continue
                key = "%s %s" % (min(u, u2), max(u, u2))
                if key not in legal:
                    legal.append(key)

        for sec in self.sections():
            values = [item[1] for item in self.items(sec)]
            for i, val in enumerate(values):
                for val2 in values[i+1:]:
                    key = "%s %s" %(min(val, val2),
                                          max(val, val2))
                    if key not in legal:
                        legal.append(key)

        for sec in self.tunnels:
            for antenna, tunnel_val in self.items(sec):
                if antenna.startswith("entrance"):
                    other_cage_antennas = self.other_cage_antenna(tunnel_val)
                    for oca in other_cage_antennas:
                        cage = self.address[oca]
                        values = [item[1] for item in self.items(sec)]
                        for cage_val in values:
                            if cage_val == tunnel_val:
                                continue
                            key = "%s %s" %(min(cage_val, tunnel_val),
                                          max(cage_val, tunnel_val))
                            if key not in legal:
                                legal.append(key)
        for l in legal:
            if l in pairs:
                pairs.remove(l)
        return pairs



        
class SetupConfig(SetupConfigMethods):
    ALL_ANTENNAS = ["1", "2", "3", "4", "5", "6", "7", "8"]
    def find_path(self, path, fname, standard, expected, possible):
        if path is None:
            self.path = data_path
            self.fname = standard
        else:
            self.path = path
            if fname is not None:
                self.fname = fname
            else:
                if os.path.isfile(os.path.join(self.path, expected)):
                    self.fname = expected
                else:
                    fnames = glob.glob(os.path.join(path, possible))
                    if len(fnames):
                        self.fname = os.path.basename(fnames[0])
                        self.path = path
                    else:
                       print("No setup config found in %s" % path)
                       self.path = data_path
                       self.fname = standard

        return os.path.join(self.path, self.fname)
        
    def __init__(self, path=None, fname=None):
        SetupConfigMethods.__init__(self)
        full_path = self.find_path(path, fname, "standard_setup.txt",
                                   'setup.txt', "setup*.txt")
        self.read(full_path)
        self.make_definitions()
        self.mismatched_pairs = self.get_mismatched_pairs()

class IdentityConfig(RawConfigParser):
    def __init__(self, path_to_fname):
        if os.path.isfile(path_to_fname):
            self.config_path = path_to_fname
        else:
            raise Exception("Could not find experiment config file %s for complex experiments",
                            path_to_fname)
        RawConfigParser.__init__(self)
        self.read(path_to_fname)

    @property
    def identity_points(self):
        field = "setup_%d_name"
        value = "point_%d_name"
        out = {}
        for section in self.sections():
            items = [item[0] for item in self.items(section)]
            setups = (len(items) - 1)//2
            if (len(items) - 1) % 2 != 0:
                raise Exception("Wrong format of ExperimentSetup config section %" % section)
            if setups < 2:
                raise Exception("Not enough setups defined in ExperimentSetup config section %" % section)

            for i in range(1, setups + 1):
                setup = self.get(section, field % i)
                point = self.get(section, value % i)
                new_key = "%s %s" % (setup, point)
                out[new_key] = self.get(section, "destination_name")
        return out


class ExperimentSetupConfig(SetupConfigMethods):
    def __init__(self, fname_with_path, **single_configs):
        """
        Read in and construct geometric description of a single experiment
        recorded using more than one EcoHAB experimental methods.

        Agrs: 
        fname_with_path: string or IdentityConfig object
           Path to the experimental setup config file, which specifies points
           (cages/tunnels) that are shared by both experimental setup


        single_configs: a dictionary of loaded SetupConfigs for each
           EcoHAB setup

        An example of an experimental setup config file
        (saved in "experiment_setup.txt"):
        [shared point 1]
        setup_1_name = ecohab_1
        point_1_name = cage A
        setup_2_name = ecohab_2
        point_2_name = cage D
        destination_name = cage A
        
        Load experiment setup config from "experiment_setup.txt", with
        config1 and config2 -- SetupConfig objects for ecohab_1 and ecohab_2:
        config = ExperimentSetupConfig(fname_with_path="experiment_setup.txt",
                                       ecohab_1=config1, ecohab_2=config2)

        """
        SetupConfigMethods.__init__(self)
        if isinstance(fname_with_path, IdentityConfig):
            experiment_config = fname_with_path
        elif isinstance(fname_with_path, str):
            experiment_config = IdentityConfig(fname_with_path)
        else:
            raise Exception("Provide a path to experiment config file or an IdentityConfig object")
        self.identity_points = experiment_config.identity_points
        self.make_sections(single_configs)
        self.ALL_ANTENNAS = self.get_all_antennas()
        self.make_definitions()

    def make_sections(self, single_configs):
        setup_names = list(single_configs.keys())
        #add individual sections from each of the setup configs
        for key in setup_names:
            this_config = single_configs[key]
            this_config_sectionss = this_config.sections()

            for section in this_config_sectionss:

                new_section_name = "%s %s" % (key, section)
                if new_section_name in self.identity_points:
                    new_section_name = self.identity_points[new_section_name]
                try:
                    self.add_section(new_section_name)
                    section_items = []
                except DuplicateSectionError:
                    section_items = [item[0] for item in self.items(new_section_name)]

                for antenna_type, value in this_config.items(section):
                    new_value = "%s_%s" %(value, key)
                    if antenna_type in section_items:
                        starts_with = antenna_type.split("_")[0]
                        how_many = len([item for item in section_items if item.startswith(starts_with)])
                        new_antenna_type = "%s_antenna%d" %(starts_with, how_many+1)
                    else:
                        new_antenna_type = antenna_type

                    self.set(new_section_name, new_antenna_type, new_value)
                    section_items = [item[0] for item in self.items(new_section_name)]
