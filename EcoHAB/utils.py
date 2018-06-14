import os

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
