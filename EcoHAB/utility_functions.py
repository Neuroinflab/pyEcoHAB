import os
from datetime import datetime
import pytz
import sys
import time

EPOCH = datetime(1970,1,1)
#UTC_OFFSET = time.mktime(EPOCH.timetuple())

EPOCH_UTC = datetime(1970, 1, 1, tzinfo=pytz.UTC)

def how_many_days(dates1, factor1):
    max_days = max([len(dates1[key]) for key in dates1.keys()])
    return max_days//factor1


def process_line_6(elements, datehourobj):
    """remove point from 2nd column of new data files"""
 
    return [elements[0],' '.join([elements[1].replace('.', ''), elements[2]]), elements[3], elements[4], elements[5]]


def process_line_7(elements, datehourobj):
    """remove point from 2nd column of new data files"""
    return [elements[0],' '.join([elements[1].replace('.', ''), elements[2]]), elements[3], elements[4], elements[5], elements[6]]


def process_line_5(elements, dateobj):
    """Add date to data (old data files)"""
    
    hour, date, datenext = dateobj
    
    if hour == '23' and elements[1][:2] == '00':
        elements[1] = ' '.join([datenetxt, elements[1]])
    else:
        elements[1] = ' '.join([date, elements[1]])
    return elements


 
def convert_time(s): 
    """Convert date and time to seconds since epoch"""
    actual_date, millisec = s.split('.')
    sec_to_epoch = time.mktime(time.strptime(actual_date, '%Y%m%d %H:%M:%S'))
    return sec_to_epoch + float(millisec)/1000

def to_timestamp_UTC(x):
    return (x - EPOCH_UTC).total_seconds()


def to_timestamp(x):
    return (x - EPOCH).total_seconds()


def is_string(obj):
    if sys.version_info < (3, 0):
        return isinstance(obj, basestring)
    else:
        return isinstance(obj, str)


def check_directory(directory, subdirectory):
    new_path = os.path.join(directory, subdirectory)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    return new_path


def results_path(path):
    head, tail = os.path.split(path)
    if tail == '':
        head, tail = os.path.split(head)
    head, tail2 = os.path.split(head)
    head = os.path.join(head, 'Results_'+tail2)
    return os.path.join(head, tail, 'Results')


def make_figure(title):
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, aspect='equal')
    fig.suptitle('%s'%(title), fontsize=14, fontweight='bold')
    return fig, ax


def make_prefix(path):
    """
    Read-in info.txt and make a prefix for all files for results.
    Parameters
    ----------
    path: str
    """
    key_list = [
        'genotype',
        'sex',
        'type of experiment',
        'date of experiment',
        'social odor',
        'no social odor',
    ]

    fname = os.path.join(path, 'info.txt')

    try:
        f = open(fname)
    except IOError:
        return ''

    info_dict = {}
    for line in f:
        try:
            key, info = line.split(':')
        except ValueError:
            continue
        info = info.strip()
        info = info.replace(' ', '_')
        info_dict[key] = info
    prefix = ''
    for key in key_list:
        if key not in info_dict:
            continue
        if key == 'social odor' or key == 'non-social odor':
            if info_dict[key] == 'none' or info_dict[key] == 'None':
                key = key.replace(' ', '_')
                prefix += '_no_' + key.replace(' ', '_') + '_'
            else:
                prefix += key.replace(' ', '_') + '_' + info_dict[key] + '_'
        else:
            prefix += info_dict[key] + '_'
    return prefix

def list_of_pairs(mice):
    pair_labels = []
    for j, mouse in enumerate(mice):
        for k in range(j+1, len(mice)):
            pair_labels.append(mice[j]+'|'+mice[k])
    return pair_labels
