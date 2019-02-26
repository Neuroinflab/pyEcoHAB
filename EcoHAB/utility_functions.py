import os
import numpy as np

def check_directory(directory,subdirectory):
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


def filter_dark(phases):
    out = []
    for phase in phases:
        if phase.endswith('dark'):
            out.append(phase)
        elif phase.endswith('DARK'):
            out.append(phase)
        elif phase.endswith('Dark'):
            out.append(phase)
    return out


def filter_light(phases):
    out = []
    for phase in phases:
        if phase.endswith('light'):
            out.append(phase)
        elif phase.endswith('LIGHT'):
            out.append(phase)
        elif phase.endswith('Light'):
            out.append(phase)
    return out


def filter_dark_light(phases):
    out = []
    for phase in phases:
        if phase.endswith('dark'):
            out.append(phase)
        elif phase.endswith('DARK'):
            out.append(phase)
        elif phase.endswith('Dark'):
            out.append(phase)
        elif phase.endswith('light'):
            out.append(phase)
        elif phase.endswith('LIGHT'):
            out.append(phase)
        elif phase.endswith('Light'):
            out.append(phase)

    return out


def add_info_mice_filename(remove_mouse):
    add_info_mice = ''
    if isinstance(remove_mouse, list):
        add_info_mice = 'remove'
        for mouse in remove_mouse:
            add_info_mice += '_' + mouse 
    elif isinstance(remove_mouse, str):
        add_info_mice = 'remove_%s' % remove_mouse
    return add_info_mice


def in_tube(antenna, next_antenna):
    if antenna % 2:
        if next_antenna  == antenna  + 1:
            return True
    else:
        if next_antenna == antenna - 1:
            return True
    return False


def in_chamber(antenna, next_antenna):
    antenna = antenna % 8
    next_antenna = next_antenna % 8
    if antenna % 2:
        if next_antenna  == antenna - 1:
            return True
    else:
        if next_antenna == antenna + 1:
            return True
    return False


def change_state(antennas):
    return np.where(abs(np.array(antennas[:-1]) - np.array(antennas[1:])) !=0)[0]

if __name__ == '__main__':
    assert in_tube(1, 2) == True
    assert in_tube(2, 1) == True
    assert in_tube(3, 2) == False
    assert in_tube(2, 3) == False
    assert in_tube(1, 8) == False
    assert in_tube(8, 1) == False
    assert in_tube(7, 8) == True
    assert in_tube(8, 7) == True
    assert in_chamber(1, 2) == False
    assert in_chamber(2, 1) == False
    assert in_chamber(3, 2) == True
    assert in_chamber(2, 3) == True
    assert in_chamber(1, 8) == True
    assert in_chamber(8, 1) == True
    assert in_chamber(7, 8) == False
    assert in_chamber(8, 7) == False
    
