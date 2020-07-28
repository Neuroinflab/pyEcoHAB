pyEcoHAB package
=====================

pyEcoHAB is a Python library for loading and analyzing behavioral data obtained in experiments performed with the Eco-HAB. Eco-HAB is a system for fully automated and ecologically relevant assessment of social impairments in mice.

pyEcoHAB supports both Python2 and Python3.

pyEcoHAB allows for easy access to antenna readings. It implements a heuristic calculating timing and duration of animal visists to EcoHAB compartments. It also provides methods for determining mouse activity (number of visits and time spent in each chamber), in-cohort sociability, solitude (time mouse spent alone in EcoHAB chambers) and following.

Installation
====================
1. Download pyEcoHAB by clicking on the green "Code" button. If you download a zip file, remember to unzip it. We add new functionality to pyEcoHAB frequently, so it is probably better to clone the repository and fetch changes every couple of weeks.
2. Add path to pyEcoHAB to your PYTHONPATH. If you are using pyCharm or Spyder, you can add it via its graphic interphase. 


Analyze data
===================
To read-in and analyze a sample data-set type:

>>> import pyEcoHAB
>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)
>>> config = pyEcoHAB.ExperimentConfigFile(pyEcoHAB.sample_data)
>>> pyEcoHAB.get_activity(data, config, 3600)
>>> pyEcoHAB.get_incohort_sociability(data, config, 3600)
>>> pyEcoHAB.get_solitude(data, config)
>>> pyEcoHAB.get_dynamic_interactions(data, config, 1000)

The standart data analysis script used by the Laboratory of Neurobiology with pyEcoHAB.sample_data as an example of the Eco-HAB dataset:

>>> import os
>>> import pyEcoHAB
>>> binsizes = [1800, 3600, 1.5*3600, 7200, 14400, 43200]
>>> path = pyEcoHAB.sample_data
>>> res_dir = os.path.join(path, "results")
>>> ehd = pyEcoHAB.Loader(path, res_dir=res_dir)
>>> config = pyEcoHAB.ExperimentConfigFile(path)
>>> pyEcoHAB.get_solitude(ehd, config)
>>> pyEcoHAB.get_incohort_sociability(ehd, config, binsize="ALL")
>>> pyEcoHAB.get_incohort_sociability(ehd, config, binsize="dark")
>>> pyEcoHAB.get_incohort_sociability(ehd, config, binsize="light")
>>> for binsize in binsizes:
...     pyEcoHAB.get_activity(ehd, config, binsize)
...     pyEcoHAB.get_incohort_sociability(ehd, config, binsize)
... 
>>> pyEcoHAB.get_dynamic_interactions(ehd, config, 1000)

Results of the data analysis can be found in pyEcoHAB/data/BALB_VPA_data_cohort_1/results/


This library is available under `GPL3 license
<http://www.gnu.org/licenses/gpl-3.0>`_.

Authors
-------
* Joanna Jędrzejewska-Szmek
* Jan Mąka
* Szymon Łęski


Acknowledgements
----------------
This software was supported by the Polish National Science Centre grant 2017/27/B/NZ4/02025.

Prerequisites
-------------
numpy and matplotlib



