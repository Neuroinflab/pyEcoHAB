import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import locale
from ExperimentConfigFile import ExperimentConfigFile
import pickle
import os

### How much time mice spend with each other

datarange = slice(0, 4, None)
datasets = ['C57 CTRL EH (3) 25.07.14', 
            'C57 CTRL EH (4) 30.07.14',
            'C57 CTRL EH (5) 20.11.14',
            'FX KO males EH (2) 19.09.14', 
            'FX WT males EH (2) 07.01.15', # 23
            ]

def splify(x):
    """Simplify a matrix: flatten and remove zeroes"""
    return x[x > 0]

for el in next(os.walk('PreprocessedData/'))[2]:
    with open('PreprocessedData/'+el, "rb") as input_file:
        ehd = pickle.load(input_file)
        ehs = pickle.load(input_file)
        path = el[:-4]
        path1 = 'RawData/' + path
        cf = ExperimentConfigFile(path1)
    phases = filter(lambda x: x.endswith('dark'), cf.sections())
    
    results = {}
    results_exp = {}
    delta = {}
    for sec in phases:
        results[sec] = splify(np.loadtxt('results_%s_%s.csv' %(path, sec), delimiter=';'))
        results_exp[sec] = splify(np.loadtxt('results_exp_%s_%s.csv' %(path, sec), delimiter=';'))
        delta[sec] = results[sec] - results_exp[sec]
    plt.figure()
    for nn in range(len(results['EMPTY 1 dark'])):
        try:
            data = map(lambda x: delta[x][nn], phases)
            idx = np.abs(np.diff(data)).sum()
        except:
            continue
        # if idx > 0.1:
        #     alpha = 0.1
        # elif idx > 0.05:
        #     alpha = 0.2
        # else:
        #     alpha = 0.8
        alpha = 0.5
        plt.plot(range(len(phases)), data, 'ko-', alpha=alpha)
        plt.gca().set_xticks(range(len(phases)))
        plt.gca().set_xticklabels(phases)
    plt.title(path)
    plt.savefig('Results/'+path+'_incohort.png')
    print '\n', path
    for key in delta:
        try:
            print key, sum(delta[key])/len(delta[key])
        except:
            print "error"

