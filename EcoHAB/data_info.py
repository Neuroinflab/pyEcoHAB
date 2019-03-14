#-*- coding: utf-8 -*-
from __future__ import print_function, division
datasets = [
    'EcoHAB_data_November/Social structure males 02.03/',
    'EcoHAB_data_November/social_dominance_swiss_webster_dominant_remove_12.02.18',
    'EcoHAB_data_November/social_structure_16.01',
    'EcoHAB_data_November/social_structure_19.01.18_rep_II',
    'EcoHAB_data_November/social_structure_swiss_webster_ctrl_05.02.18',
    'EcoHAB_data_November/mice K Wisniewska',
    'EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/',
    'EcoHAB_data_November/C57 13-24.04 long/',
    'EcoHAB_data_November/C57 males long 11-22.05.18/',
    'EcoHAB_data_November/C57 TIMP rep 2/',
    'EcoHAB_data_November/C57 males rep 2/',
    'EcoHAB_data_November/C57 males TIMP/',
    'EcoHAB_data_November/BTBR males/',
    'EcoHAB_data_November/long_experiment_WT',
    'EcoHAB_data_November/Maciek_01_30_2018',
    'EcoHAB_data_November/Eksperyment_A.Krzemień_22-25.05.2018/',
    "EcoHAB_data_November/C57_males_long_26.05-06.06.2018_after_TIMP/",
    "EcoHAB_data_November/mice_K_Wisniewska_short_familiar/",
    "EcoHAB_data_November/FXKO_F_S_EHcl_familiar/",
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18",
    "EcoHAB_data_November/FXKO_F_10_L_EHcl_12.06-22.06_rep1",
    "EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark/",
    "EcoHAB_data_November/long_experiment_KO_from_phase_SNIFF_10_dark/",
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/",
    "EcoHAB_data_November/WT_long_03-13.07.2018_Asia/",
    "EcoHAB_data_November/males_social_dominance_21-25.07.2018-Maciek/",
    "EcoHAB_data_November/Asia_krzemień_kontrola_27.07-01.08.2018/",
    "EcoHAB_data_November/C57 social contagion stress/",
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18",
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/",
    "EcoHAB_data_November/C57 social contagion stress/",
    "EcoHAB_data_November/Asia_WT_FX_mixed_13-24.08.2018/",
    "EcoHAB_data_November/Asia_Krzemień_AL1KO_samice_08.10-11.10/",
    "EcoHAB_data_November/Asia_Krzemiń_Kontrola_samice_02.10-06.10",
    "EcoHAB_data_November/Asia_Krzemień_Kontrola_samce_II_28.09-01.10",
    "EcoHAB_data_November/C57_males_social interaction-12.10-22.10/",
    "EcoHAB_data_November/C57_males_emotional_contagion-22.10-26.10/",
    "EcoHAB_data_November/C57_males_reward_5-9.11",
    "EcoHAB_data_November/C57_males_emotional_contagion_stress-12.11-16.11",
    "EcoHAB_data_November/Myszy_zlote_kizinski_03-07.12",
    "EcoHAB_data_November/Koziński_turkusowe_10-14.12/",
    "EcoHAB_data_November/C57_males_before_TIMP_Social_Interaction_21.12-31.12.18/",
    "EcoHAB_data_November/C57_males_emotional_contagion_reward_31.12-04.01",
    "EcoHAB_data_November/Asia_Krzemień_Al1KO_04-08.01.19/",
    "EcoHAB_data_November/C57_long_afterTIMP_08-17.01.19/",
    "EcoHAB_data_November/Turkusowe II 04-08.01.2019/",
    "EcoHAB_data_November/Myszy z Krakowa/",
    "EcoHAB_data_November/BTBR_LONG_KSENIA_12-22.02.19/",
    "EcoHAB_data_November/Samce C57 sham 3 rep 22.02-04.03/",
    "EcoHAB_data_November/Samce C57 sham EH dominacja 05-08.03.2019/",
    "EcoHAB_data_November/kraków dominacja 08.03-12.03.2019/"
]

remove_tags = {
    'EcoHAB_data_November/C57 males rep 2/':['0065-0161984735'],
    'EcoHAB_data_November/BTBR males/':['0065-0136658439',
                                                     '0065-0141855614'],
    # "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18":    ['0065-0136673085', '0065-0136660668'],
    # "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/":    ['0065-0136673085', '0065-0136660668'],
    # "EcoHAB_data_November/C57 social contagion stress/":    ['0065-0136673085', '0065-0136660668']
}
how_many_appearances = {
    'EcoHAB_data_November/C57 males rep 2/':1000,
    'EcoHAB_data_November/BTBR males/':500,
    'EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/':200,
    "EcoHAB_data_November/males_social_dominance_21-25.07.2018-Maciek/":200,
    "EcoHAB_data_November/C57 males reward-neutral 26.06-01.07.18": 200,
    "EcoHAB_data_November/C57 males reward-neutral control 16.07-20.07.18/": 200,
    "EcoHAB_data_November/C57 social contagion stress/": 200,
}
antenna_positions = {
    'EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark/':{'1': 1,
                                                                                          '2': 5,
                                                                                          '3': 3,
                                                                                          '4': 6,
                                                                                          '5': 4,
                                                                                          '6': 2,
                                                                                          '7': 7,
                                                                                          '8': 8}}
binsizes = [12 * 3600., 2 * 3600.,1*3600.]
bintitles = ['12', '2']
standard_ant_pos = {'1': 5, '2': 6, '3': 7, '4': 8, '5': 1, '6': 2, '7': 3, '8': 4}
basic = ['Number of visits to box %d\n','Total time in box %d, seconds\n']
standard_cages = {'1': 1, '2': 2, '3': 3, '4': 4}
non_standard_cages = {}
standard_headers = {}
for i in range(1,5):
    standard_headers[str(i)] = basic
non_standard_headers = {}

