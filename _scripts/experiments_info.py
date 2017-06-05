#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:11:02 2016

@author: jmaka
"""

datasets = ['BALB EH (1) 06.03.14', # 0
            'BALB EH (2) po raz 2 10.04.14',
            'C57 EH (1) 23.04.14',
            'FX KO F (3) EH (1) 15.05.14',
            'BALB VPA EH (1) 05.05.14', # ***
            'BALB VPA EH (2) 09.05.14', # 5
            'C57 CTRL EH (2) 28.04.14',
            'BTBR 1st time - stary EH',
            'EH BTBR po raz 2', #8
            'BALB CTRL EH (3) 12.06.14 - 15.06.14',
            'BALB VPA EH (4) 16.06.14 - 19.06.14', # 10
            'C57 VPA EH (1) 20.06.14',
            'C57 CTRL EH (3) 25.07.14', # 12
            'C57 CTRL EH (4) 30.07.14',
            'FX KO males (1) 08.08.14',
            'C57 VPA EH (1Biol) 28.10.2014', # 15, nowy Eco-HAB
            'BALB CTRL EH (4) 09.09.14',
            'BALB CTRL EH (5) 15.09.14',
            'FX WT males EH (1) 24.09.14',
            'FX KO females EH 2 29.09.14',
            'FX KO males EH (2) 19.09.14', # 20
            'C57 CTRL EH (5) 20.11.14',
            'FX KO males EH (3) 14.01.15', # 22
            'FX WT males EH (2) 07.01.15', # 23
            'FX WT females EH (1) 16.02.15', # 24
            'FX WT females EH (2) 20.02.15',
            'FX KO females EH (3) 11.02.15',
            'BTBR females (2) EH - 13.03.2015',
            'BALB VPA EH Biol (1) 23.03.15',
            'C57 VPA EH Biol (2) 27.03.15',
            'BTBR males EH (2) 31.03.15', # 30
            'C57 NaCl Biol (1) 26.06.15',
            'C57 NaCl Biol (2) 01.07.15',
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1',
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2',
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1',
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2',
            'BALB NaCl (1) Rep (1) 09.09.15',
            'BALB NaCl (2) Rep (1) 02.10.15',
            'BALB NaCl (2) Rep (2) 11.10.15', # 39
            'FX WT females EH Nen (1) 07.12.15', 
            'FAMB KO males (1) rep 2 26.02.16',
            'FAMB WT males (1) rep 1 22.02.16',
            'FAMB WT males (1) rep 2 02.03.16',
            'FX WT females EH (B) 05.05.16', # 44
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16',
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16',
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16',
            'PV Cre males HT EH (1) 09.06.16', # 50
            'test_KO']
for i in datasets:
    print i
smells = {'BALB EH (1) 06.03.14': {'soc': 3, 'nsoc': 1},
            'BALB EH (2) po raz 2 10.04.14': {'soc': 3, 'nsoc': 1},
            'C57 EH (1) 23.04.14': {'soc': 3, 'nsoc': 1},
            'FX KO F (3) EH (1) 15.05.14': {'soc': 3, 'nsoc': 1},
            'BALB VPA EH (1) 05.05.14': {'soc': 4 , 'nsoc': 2}, 
            'BALB VPA EH (2) 09.05.14': {'soc': 4, 'nsoc': 2},
            'C57 CTRL EH (2) 28.04.14': {'soc': 4, 'nsoc': 2},
            'BTBR 1st time - stary EH': {'soc': 3, 'nsoc': 1},
            'EH BTBR po raz 2': {'soc': 3, 'nsoc': 1},
            'BALB CTRL EH (3) 12.06.14 - 15.06.14': {'soc': 1, 'nsoc': 2}, # soc miedzy 2, 3; nsoc 4, 5
            'BALB VPA EH (4) 16.06.14 - 19.06.14': {'soc': 3, 'nsoc': 1}, # soc miedzy 6, 7; nsoc 2, 3
            'C57 VPA EH (1) 20.06.14': {'soc': 3, 'nsoc': 4},
            'C57 CTRL EH (3) 25.07.14': {'soc': 1, 'nsoc': 3},
            'C57 CTRL EH (4) 30.07.14': {'soc': 3, 'nsoc': 1},
            'FX KO males (1) 08.08.14': {'soc': 1, 'nsoc': 3},
            'C57 VPA EH (1Biol) 28.10.2014': {'soc': 1, 'nsoc': 3},
            'BALB CTRL EH (4) 09.09.14': {'soc': 1, 'nsoc': 3},
            'BALB CTRL EH (5) 15.09.14': {'soc': 3, 'nsoc': 1},
            'FX WT males EH (1) 24.09.14': {'soc': 1, 'nsoc': 3},
            'FX KO females EH 2 29.09.14': {'soc': 3, 'nsoc': 1},
            'FX KO males EH (2) 19.09.14': {'soc': 3, 'nsoc': 1},
            'C57 CTRL EH (5) 20.11.14': {'soc': 3, 'nsoc': 1},
            'FX KO males EH (3) 14.01.15': {'soc': 1, 'nsoc': 2},
            'FX WT males EH (2) 07.01.15': {'soc': 1, 'nsoc': 2}, # 23
            'FX WT females EH (1) 16.02.15': {'soc': 3, 'nsoc': 1}, # 24
            'FX WT females EH (2) 20.02.15': {'soc': 1, 'nsoc': 3},
            'FX KO females EH (3) 11.02.15': {'soc': 1, 'nsoc': 3},
            'BTBR females (2) EH - 13.03.2015': {'soc': 1, 'nsoc': 3},
            'BALB VPA EH Biol (1) 23.03.15': {'soc': 3, 'nsoc': 1},
            'C57 VPA EH Biol (2) 27.03.15': {'soc': 3, 'nsoc': 1},
            'BTBR males EH (2) 31.03.15': {'soc': 3, 'nsoc': 1},
            'C57 NaCl Biol (1) 26.06.15': {'soc': 1, 'nsoc': 3},
            'C57 NaCl Biol (2) 01.07.15': {'soc': 3, 'nsoc': 1},
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1': {'soc': 3, 'nsoc': 1},
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2': {'soc': 1, 'nsoc': 3},
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1': {'soc': 1, 'nsoc': 3},
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2': {'soc': 3, 'nsoc': 1},
            'BALB NaCl (1) Rep (1) 09.09.15': {'soc': 3, 'nsoc': 1},
            'BALB NaCl (2) Rep (1) 02.10.15': {'soc': 1, 'nsoc': 3},
            'BALB NaCl (2) Rep (2) 11.10.15': {'soc': 1, 'nsoc': 3},
            'FX WT females EH Nen (1) 07.12.15': {'soc': 1, 'nsoc': 3},
            'FAMB KO males (1) rep 2 26.02.16': {'soc': 1, 'nsoc': 3},
            'FAMB WT males (1) rep 1 22.02.16': {'soc': 1, 'nsoc': 3},
            'FAMB WT males (1) rep 2 02.03.16': {'soc': 3, 'nsoc': 1},
            'FX WT females EH (B) 05.05.16': {'soc': 3, 'nsoc': 1},
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16': {'soc': 3, 'nsoc': 1},
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16': {'soc': 1, 'nsoc': 3},
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16': {'soc': 3, 'nsoc': 1},
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16': {'soc': 3, 'nsoc': 1},
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16': {'soc': 1, 'nsoc': 3},
            'PV Cre males HT EH (1) 09.06.16': {'soc': 3, 'nsoc': 1},
            'Random_test':{'soc': 1, 'nsoc': 3},
            }
antenna_positions = {'BALB EH (1) 06.03.14': None,
            'BALB EH (2) po raz 2 10.04.14': None,
            'C57 EH (1) 23.04.14': {'6': 1, '5': 2, '7': 3, '8': 4, '2': 5, '1': 6, '4': 7, '3': 8},
            'FX KO F (3) EH (1) 15.05.14': None,
            'BALB VPA EH (1) 05.05.14': None,
            'BALB VPA EH (2) 09.05.14': None,
            'C57 CTRL EH (2) 28.04.14': None, 
            'BTBR 1st time - stary EH': None,
            'EH BTBR po raz 2': None,
            'BALB CTRL EH (3) 12.06.14 - 15.06.14': None,
            'BALB VPA EH (4) 16.06.14 - 19.06.14': None,
            'C57 VPA EH (1) 20.06.14': None,
            'C57 CTRL EH (3) 25.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
            'C57 CTRL EH (4) 30.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
            'FX KO males (1) 08.08.14': {'5': 1, '6': 2, '7': 3, '8': 4, '1': 5, '2': 6, '3': 7, '4': 8},
            'C57 VPA EH (1Biol) 28.10.2014': None,
            'BALB CTRL EH (4) 09.09.14': None,
            'BALB CTRL EH (5) 15.09.14': None,
            'FX WT males EH (1) 24.09.14': None,
            'FX KO females EH 2 29.09.14': None,
            'FX KO males EH (2) 19.09.14': None,
            'C57 CTRL EH (5) 20.11.14': None,
            'FX KO males EH (3) 14.01.15': None,
            'FX WT males EH (2) 07.01.15': None,
            'FX WT females EH (1) 16.02.15': None, # 24
            'FX WT females EH (2) 20.02.15': None,
            'FX KO females EH (3) 11.02.15': None,
            'BTBR females (2) EH - 13.03.2015': None,
            'BALB VPA EH Biol (1) 23.03.15': None,
            'C57 VPA EH Biol (2) 27.03.15': None,
            'BTBR males EH (2) 31.03.15': None,
            'C57 NaCl Biol (1) 26.06.15': None,
            'C57 NaCl Biol (2) 01.07.15': None,
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1': None,
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2': None,
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1': None,
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2': None,
            'BALB NaCl (1) Rep (1) 09.09.15': None,
            'BALB NaCl (2) Rep (1) 02.10.15': None,
            'BALB NaCl (2) Rep (2) 11.10.15': None,
            'FX WT females EH Nen (1) 07.12.15': None,
            'FAMB KO males (1) rep 2 26.02.16': None,
            'FAMB WT males (1) rep 1 22.02.16': None,
            'FAMB WT males (1) rep 2 02.03.16': None,
            'FX WT females EH (B) 05.05.16': None,
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16': None,
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16': None,
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16': None,
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16': None,
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16': None,
            'PV Cre males HT EH (1) 09.06.16': None,
            'Random_test': None,
            }