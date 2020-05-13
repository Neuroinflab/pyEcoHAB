#!/usr/bin/env python
# encoding: utf-8
from __future__ import division, absolute_import, print_function
"""
ExperimentConfigFile.py

Created by Szymon Łęski on 2013-02-19.
Copyright (c) 2013 Laboratory of Neuroinformatics. All rights reserved.
"""

import os    
import numpy as np                                           
import sys
import time
if sys.version_info < (3, 0):
    from ConfigParser import RawConfigParser, NoSectionError
else:
    from configparser import RawConfigParser, NoSectionError

import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')

import matplotlib.ticker
import matplotlib.dates as mpd
import matplotlib.pyplot as plt



class ExperimentConfigFile(RawConfigParser, matplotlib.ticker.Formatter):
    def __init__(self, path, fname=None):    
        RawConfigParser.__init__(self)
        self.path = path               
        if fname is None:
            if os.path.isfile(os.path.join(path, 'config.txt')):
                self.fname = 'config.txt'
            else:
                self.fname = filter(lambda x: x.startswith('config') 
                        and x.endswith('.txt'), os.listdir(path))[0]
        else:                  
            self.fname = fname
        self.read(os.path.join(path, self.fname)) 
   
    def gettime(self, sec): 
        """Convert start and end time and date read from section sec
        (might be a list)
        of the config file 
        to a tuple of times from epoch."""
        if type(sec) == list:
            starts = []
            ends = []
            for ss in sec:
                st, et = self.gettime(ss)
                starts.append(st)
                ends.append(et)
            return min(starts), max(ends)
        else:
            tstr1 = self.get(sec, 'startdate') + self.get(sec, 'starttime')
            tstr2 = self.get(sec, 'enddate') + self.get(sec, 'endtime')
            if len(tstr1) == 15:
                t1 = time.strptime(tstr1, '%d.%m.%Y%H:%M')
            elif len(tstr1) == 18:                        
                t1 = time.strptime(tstr1, '%d.%m.%Y%H:%M:%S')
            else: 
                raise Exception('Wrong date format in %s' %self.fname)

            if len(tstr2) == 15:
                t2 = time.strptime(tstr2, '%d.%m.%Y%H:%M')
            elif len(tstr2) == 18:                        
                t2 = time.strptime(tstr2, '%d.%m.%Y%H:%M:%S')
            else: 
                raise Exception('Wrong date format in %s' %self.fname)

            return time.mktime(t1), time.mktime(t2)

    def __call__(self, x, pos=0):
        x = mpd.num2epoch(x)
        for sec in self.sections():
            t1, t2 = self.gettime(sec)
            if t1 <= x and x < t2:
                return sec
        return 'Unknown'    
    
    def mark(self, sec, ax=None):
        """Mark given phases on the plot"""
        if ax is None:
            ax = plt.gca()
        ylims = ax.get_ylim()
        for tt in self.gettime(sec):
            ax.plot([mpd.epoch2num(tt),] * 2, ylims, 'k:')
        plt.draw()
    
    def plot_nights(self, sections, ax=None):
        """Plot night from sections"""
        if ax is None:
            ax = plt.gca()
        ylims = ax.get_ylim()  
        xlims = ax.get_xlim()
        if type(sections) == str:
            sections = [sections]  
        for sec in sections:
            t1, t2 = self.gettime(sec)        
            plt.bar(mpd.epoch2num(t1), ylims[1] - ylims[0], 
                    width=mpd.epoch2num(t2) - mpd.epoch2num(t1), 
                    bottom=ylims[0], color='0.8', alpha=0.5, zorder=-10)
        ax.set_xlim(xlims)
        plt.draw()
    
    def plot_sections(self):
        """Diagnostic plot of sections defined in the config file."""
        figg = plt.figure()                         
        for idx, sec in enumerate(self.sections()):
            t1, t2 = mpd.epoch2num(self.gettime(sec)) #cf2time(cf, sec)
            plt.plot([t1, t2], [idx, idx], 'ko-') 
            plt.plot([t2], [idx], 'bo')
            plt.text(t2 + 0.5, idx, sec)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpd.HourLocator(np.array([00]), 
                                                   tz=tzone)) 
        ax.xaxis.set_major_formatter(mpd.DateFormatter('%d.%m %H:%M', tz=tzone))
        ax.autoscale_view()
        ax.get_figure().autofmt_xdate()
        plt.title(self.path) 
        plt.draw()
        
