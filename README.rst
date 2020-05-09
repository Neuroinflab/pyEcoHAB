pyEcoHAB package
=====================

pyEcoHAB is a Python library for loading and analyzing behavioral data obtained in experiments performed with the Eco-HAB. Eco-HAB is a system for fully automated and ecologically relevant assessment of social impairments in mice.

pyEcoHAB supports both Python2 and Python3.

pyEcoHAB allows for easy access to antenna readings. It implements a heuristic calculating timing and duration of animal visists to EcoHAB compartments. It also provides methods for determining mouse activity (number of visits and time spent in each chamber), in-cohort sociability, solitude (time mouse spent alone in EcoHAB chambers) and following.

To read-in and analyze a sample data-set type:
>>> import pyEcoHAB
>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)
>>> config = pyEcoHAB.ExperimentConfigFile(pyEcoHAB.sample_data)
>>> pyEcoHAB.get_activity(data, config, 3600)
>>> pyEcoHAB.get_incohort_sociability(data, config, 3600)
>>> pyEcoHAB.get_solitude(data, config)
>>> pyEcoHAB.get_following(data, config, 1000)


The library is available under `GPL3 license
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



