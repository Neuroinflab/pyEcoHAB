 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 14:56:11 2017

@author: jmaka
"""


rawdatasets = [
            'BALB EH (1) 06.03.14,
            'BALB EH (2) po raz 2 10.04.14',
            'C57 EH (1) 23.04.14',
            'FX KO F (3) EH (1) 15.05.14',
            'BALB VPA EH (1) 05.05.14', # ***
            'BALB VPA EH (2) 09.05.14', # 5
            'C57 CTRL EH (2) 28.04.14',
            'BTBR 1st time - stary EH',
            'EH BTBR po raz 2', #8
            'BALB CTRL EH (3) 12.06.14 - 15.06.14',
            'BALB VPA EH (4) 16.06.14 - 19.06.14', # 10 #
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
            ]
            
dataset = []





#{"id":0,"fenotype":"BALB", "type":"EH","date"="06.03.14", "group_nr":1, "comment":""}
paths_WT = [#'FX WT females EH Nen (1) 07.12.15',
            'FX WT females EH (1) 16.02.15',
            'FX WT females EH (B) 05.05.16',
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16',
            #'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16',
            ]
paths_KO =  ['FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16',
            #'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16',
            #'FX KO F (3) EH (1) 15.05.14',
            'FX KO females EH 2 29.09.14',
            'FX KO females EH (3) 11.02.15'
            ]
             
paths_ = [              
          'C57 VPA EH (1) 20.06.14'
          'C57 VPA EH (1Biol) 28.10.2014', # 15, nowy Eco-HAB
          'C57 VPA EH Biol (2) 27.03.15',
         ]
         

paths_ = [  'C57 CTRL EH (3) 25.07.14', # 12
            'C57 CTRL EH (4) 30.07.14',
            'C57 CTRL EH (2) 28.04.14',
         ]      
         
paths_ = [
          'BALB EH (1) 06.03.14', # 0
          'BALB EH (2) po raz 2 10.04.14',
          
          ]
          
          

