from __future__ import division, print_function
import os
import analiza_friends as af
import cage_visits as cv
import mouse_speed as ms
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile


binsizes = [12 * 3600., 2 * 3600.]
bintitles = ['12', '2']
standard_ant_pos = {'1': 5, '2': 6, '3': 7, '4': 8, '5': 1, '6': 2, '7': 3, '8': 4}
basic = ['Number of visits to box %d\n','Total time in box %d, seconds\n']
standard_cages = {'1': 1, '2': 2, '3': 3, '4': 4}
non_standard_cages = {}
standard_headers = {}
for i in range(1,5):
    standard_headers[str(i)] = basic
non_standard_headers = {}

datasets = [
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18",
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/",
    "EcoHAB_data_November/C57 social contagion stress/"
]

remove_tags = {
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18":
    ['0065-0136673085', '0065-0136660668'],
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/":
    ['0065-0136673085', '0065-0136660668'],
    "EcoHAB_data_November/C57 social contagion stress/":
    ['0065-0136673085', '0065-0136660668']
}

how_many_appearances = {
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18": 200,
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/": 200,
    "EcoHAB_data_November/C57 social contagion stress/": 200,
}
antenna_positions = {}
homepath = os.path.expanduser("~/")
threshold = 3
if __name__ == '__main__':
    

    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
            if isinstance(remove_mouse, list):
                which_mice = ''
                for rm in remove_mouse:
                    which_mice += rm+'_'
            elif isinstance(remove_mouse, str):
                which_mice = remove_mouse
        else:
            remove_mouse = None
            which_mice = ''

        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 500
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path)
        af.get_mouse_alone(ehs, cf)
        af.get_in_cohort_sociability(ehs, cf,
                                     remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, which_phases="ALL",
                                 remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, which_phases="dark",
                                 remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, which_phases="light",
                                 remove_mouse=remove_mouse)
        ms.get_following(ehs, cf)
        for binsize in binsizes:
            print('Binsize ',binsize/3600)
            cv.get_visits(ehs, cf)
