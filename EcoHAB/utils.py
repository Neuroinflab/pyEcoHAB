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
    head = os.path.join(head,'Results_'+tail2)
    return os.path.join(head,tail)

def make_figure(title):
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, aspect='equal')
    fig.suptitle('%s'%(title), fontsize=14, fontweight='bold')
    return fig, ax
    
