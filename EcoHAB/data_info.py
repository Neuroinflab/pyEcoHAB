# -*- coding: utf-8 -*-
datasets = [
    # '/home/jszmek/EcoHAB_data_November/Social structure males 02.03/',
    # '/home/jszmek/EcoHAB_data_November/social_dominance_swiss_webster_dominant_remove_12.02.18',
    # '/home/jszmek/EcoHAB_data_November/social_structure_16.01',
    # '/home/jszmek/EcoHAB_data_November/social_structure_19.01.18_rep_II',
    # '/home/jszmek/EcoHAB_data_November/social_structure_swiss_webster_ctrl_05.02.18',
    # '/home/jszmek/EcoHAB_data_November/mice K Wisniewska',
    # '/home/jszmek/EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/',
    # '/home/jszmek/EcoHAB_data_November/C57 13-24.04 long/',
    # '/home/jszmek/EcoHAB_data_November/C57 males long 11-22.05.18/',
    # '/home/jszmek/EcoHAB_data_November/C57 TIMP rep 2/',
    # '/home/jszmek/EcoHAB_data_November/C57 males rep 2/',
    # '/home/jszmek/EcoHAB_data_November/C57 males TIMP/',
    # '/home/jszmek/EcoHAB_data_November/BTBR males/',
    # '/home/jszmek/EcoHAB_data_November/long_experiment_WT',
    # '/home/jszmek/EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark',
    # '/home/jszmek/EcoHAB_data_November/long_experiment_KO_from_phase_SNIFF_10_dark',
    # '/home/jszmek/EcoHAB_data_November/Maciek_01_30_2018',
    '/home/jszmek/EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/',
    '/home/jszmek/EcoHAB_data_November/Eksperyment_A.Krzemień_22-25.05.2018/',
]

remove_tags = {
    '/home/jszmek/EcoHAB_data_November/C57 males rep 2/':['0065-0161984735'],
    '/home/jszmek/EcoHAB_data_November/BTBR males/':['0065-0136658439',
                                                     '0065-0141855614']
}
how_many_appearances = {
    '/home/jszmek/EcoHAB_data_November/C57 males rep 2/':1000,
    '/home/jszmek/EcoHAB_data_November/BTBR males/':500,
    '/home/jszmek/EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/':200
}
antenna_positions = {
    '/home/jszmek/EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark':{'1': 1,
                                                                                                       '2': 5,
                                                                                                       '3': 3,
                                                                                                       '4': 6,
                                                                                                       '5': 4,
                                                                                                       '6': 2,
                                                                                                       '7': 7,
                                                                                                       '8': 8}}
binsizes = [12 * 3600., 1 * 3600.]
bintitles = ['12', '2']
standard_ant_pos = {'1': 5, '2': 6, '3': 7, '4': 8, '5': 1, '6': 2, '7': 3, '8': 4}
basic = ['Number of visits to box %d\n','Total time in box %d, seconds\n']
standard_cages = {1: 1, 2: 2, 3: 3, 4: 4}
non_standard_cages = {}
standard_headers = {}
for i in range(1,5):
    standard_headers[i] = basic
non_standard_headers = {}

