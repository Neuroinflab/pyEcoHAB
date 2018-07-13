from __future__ import division, print_function
import os
import utils
import numpy as np

def write_single_chamber(f, header, heads, address, mice, phases, time, data_stim ):
    
    for j, h in enumerate(heads):
        f.write(h%address)
        f.write(header+'\n')
        for i, mouse in enumerate(mice):
            
            lines = [mouse for i in data_stim[j][phases[0]][mouse]]
            
            for phase in phases:
                for k,t in enumerate(time[phase]):
                    if phase == phases[0]:
                        lines[k] += ';%3.2f'%(t/3600)

                    lines[k] += ';'+str(data_stim[j][phase][mouse][k])
                                                                
            for line in lines:
                f.write(line+'\n')

def save_data_cvs(data, fname, path, which, headers):
    new_path = os.path.join(path, 'time_in_chambers')
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    fname = os.path.join(new_path, fname)
    f = open(fname,'w')
    phases = data['phases']
    mice = data['mice']
    
    header = 'mouse;\"time [h]\"'
    for phase in phases:
        header += ';\"'+phase+"\""
  
    for stim in which.keys():
        heads = headers[stim]
        write_single_chamber(f,header, heads, which[stim], mice, phases, data['time'], data[stim])


def save_single_histograms(result, fname, mice, phase, main_directory, directory, prefix, additional_info=None):
    new_name = os.path.join(directory, 'data')
    directory = utils.check_directory(main_directory, new_name)
    if isinstance(additional_info, str):
        fname =  os.path.join(directory, '%s_%s_%s_%s.csv'% (fname, prefix, phase, additional_info))
    else:
        fname =  os.path.join(directory, '%s_%s_%s.csv'% (fname, prefix, phase))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
        return None
    for i, mouse in enumerate(mice):
        f.write(';'+mouse)
    f.write(';\n')
    for i, mouse in enumerate(mice):
        f.write(mouse + ';')
        for j, mouse in enumerate(mice):
            f.write(str(result[i, j])+';')
        f.write('\n')


def write_csv_rasters(mice, phases, output, directory, dirname, fname):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(directory, new_name)
    fname = os.path.join(directory, fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
   
    header = 'mouse pair'
    for phase in phases:
        header += ';' + phase
        
    header += '\n'
    f.write(header)
    new_output, pairs = make_table_of_pairs(output, phases, mice)
    for i, pair in enumerate(pairs):
        f.write(pair)
        for j in range(len(phases)):
            f.write(';')
            f.write(str(new_output[i,j]))
        f.write('\n')
    f.close()

def write_csv_tables(results, phases, mice, main_directory, dirname, fname, prefix):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s.csv' % (fname, prefix))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    phase_pairs = [phases[2*i]+(len(mice)+2)*';'+ phases[2*i+1]+'\n' for i in range(len(phases)//2)]
    new_results = [(results[2*i], results[2*i+1])for i in range(len(phases)//2)]

    mice_header = ';'
    
    for i in range(2):
        for mouse in mice:
            mice_header += mouse + ';'
        if not i:
            mice_header += ';;'
    mice_header += '\n'
    
    if len(phases) % 2:
        phase_pairs.append((phases[-1]+'\n'))
        new_results.append((results[-1],))
    for i, new_pair in enumerate(phase_pairs):
        f.write(new_pair)
        if len(new_results[i]) > 1:
            f.write(mice_header)
        else:
            single = 1+len(mice)*(len(mice[0])+1)
            f.write(mice_header[:single]+'\n')
        for j, mouse in enumerate(mice):
            for l, nr in enumerate(new_results[i]):
                f.write(mouse+';')
                for k, mouse in enumerate(mice):
                    f.write(str(nr[j, k])+';')
                if not l:
                    f.write(';')
            f.write('\n')
        f.write('\n')
    f.close()

def write_csv_alone(alone, phases, mice, main_directory, prefix):
    directory = utils.check_directory(main_directory, 'mouse_alone')
    fname =  os.path.join(directory, 'mouse_alone_%s.csv' % prefix)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    header = 'Mice alone in chamber %d\n'
    phases_header = 'Mouse;'
    for phase in phases:
        phases_header += phase + ';'
    phases_header += '\n'
    for i in range(1, 5):
        f.write(header % i)
        f.write(phases_header)
        for j, mouse in enumerate(mice):
            f.write(mouse+';')
            for k, phase in enumerate(phases):
                f.write(str(alone[i-1, j, k])+';')
            f.write('\n')
    f.close()
    
def make_table_of_pairs(FAM, phases, mice):
    new_shape = (len(mice)*(len(mice)-1)//2, len(phases))
    output = np.zeros(new_shape)
    pair_labels = utils.list_of_pairs(mice)
    for i, phase in enumerate(phases):
        l = 0
        for j, mouse in enumerate(mice):
            for k in range(j+1, len(mice)):
                output[l,i] = FAM[i,j,k]
                l += 1

    return output, pair_labels


