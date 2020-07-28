pyEcoHAB package
=====================

pyEcoHAB is a Python library for loading and analyzing behavioral data obtained in experiments performed with the Eco-HAB. Eco-HAB is a system for fully automated and ecologically relevant assessment of social impairments in mice.

pyEcoHAB supports both Python2 and Python3.

pyEcoHAB allows for easy access to antenna readings. It implements a heuristic calculating timing and duration of animal visists to EcoHAB compartments. It also provides methods for determining mouse activity (number of visits and time spent in each chamber), in-cohort sociability, solitude (time mouse spent alone in EcoHAB chambers) and following.

Install
To read-in and analyze a sample data-set type:

>>> import pyEcoHAB
>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)
>>> config = pyEcoHAB.Timeline(pyEcoHAB.sample_data)
>>> pyEcoHAB.get_activity(data, config, 3600)
>>> pyEcoHAB.get_incohort_sociability(data, config, 3600)
>>> pyEcoHAB.get_solitude(data, config)
>>> pyEcoHAB.get_dynamic_interactions(data, config, 1000)

Standart data analysis script used by the Laboratory of Neurobiology with pyEcoHAB.sample_data as an example of an Eco-HAB dataset:

>>> import os
>>> import pyEcoHAB
>>> binsizes = [1800, 3600, 1.5*3600, 7200, 14400, 43200]
>>> path = pyEcoHAB.sample_data
>>> res_dir = os.path.join(path, "results")
>>> ehd = pyEcoHAB.Loader(path, res_dir=res_dir)
>>> config = pyEcoHAB.Timeline(path)
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

Analyze a custom Eco-HAB experiment
=========

If you want to analyse data obtained by a non-standard EcoHAB setup or a modular experimental setup utilizing more than one EcoHAB setup, you should provide a config file for every non-standard EcoHAB setup and a separate config file describing how EcoHAB setups are combined together into a modular experimental setup.

Config file for an Eco-HAB setup
------------
pyEcoHAB requires a setup.txt file describing geometry of the experimental setup and antenna positions in the data directory. If you don't provide one, pyEcoHAB will load a configuration of a standard EcoHAB setup with 4 chambers connected with 4 tunnels, forming a square with a chamber in each corner and tunnels providing sides of that square. A standard setup configuration file provided by pyEcoHAB:


   [setup]

   name = default

   
   [cage B]

   entrance_antenna1 = 2

   entrance_antenna2 = 3

   
   [cage C]

   entrance_antenna1 = 4

   entrance_antenna2 = 5

   
   [cage D]

   entrance_antenna1 = 6

   entrance_antenna2 = 7

   
   [cage A]

   entrance_antenna1 = 8

   entrance_antenna2 = 1


   [tunnel 1]

   entrance_antenna1 = 1

   entrance_antenna2 = 2


   [tunnel 2]

   entrance_antenna1 = 3

   entrance_antenna2 = 4


   [tunnel 3]

   entrance_antenna1 = 5

   entrance_antenna2 = 6


   [tunnel 4]

   entrance_antenna1 = 7

   entrance_antenna2 = 8

   
In a setup.txt file you need to specify:
a. Your setups name in section [setup]. This is very important for modular EcoHAB setups, because every submodule has to have a unique name.
b. Chambers of the setup and their entrance antennas and internal antennas (if there are any internal antennas). You add each chamber specification as a separate section (in square brackets). Every chamber name needs to be unique and contain the word cage (lower case). In the section specifying each chamber list entrance antennas and internal antennas and their numbers. If there is more than one  antenna of a certain type you need to number them e.g. external_antenna1 = 6, external_antenna2 = 7.
c. tunnels connecting chambers.  You add each tunnel specification as a separate section (in square brackets). Every tunnel name needs to be unique and contain the word tunnel (lower case). In the section specifying each tunnel list entrance antennas and internal antennas and their numbers.  If there is more than one  antenna of a certain type you need to number them e.g. external_antenna1 = 7, external_antenna2 = 8.

A configuration file for a custom setup with two chambers connected with a tunnel with an additional internal antenna in cage A: 

   [setup] 

   name = my_experiment
      

   [cage A]

   external_antenna = 1

   internal_antenna = 3


   [cage B]

   external_antenna = 2


   [tunnel 1]

   external_antenna1 = 1
   
   external_antenna2 = 2

Config file for a modular Eco-HAB setup
------------
If your experimental setup consists of more then one Eco-HAB experimental setups, you need to provide a setup config file for every setup and a master configuration setup file describing the whole setups and mainly what chambers/tunnels were parts of at least two setups.

Example 1
~~~~~~~~~
An experiment consisiting of a standard Eco-HAB setup with additional internal antennas in cage A (antenna 1) and cage C (antenna 8):
1. Standard Eco-Hab setup can be provided by pyEcoHAB:
   [setup]

   name = default

   
   [cage B]

   entrance_antenna1 = 2

   entrance_antenna2 = 3

   
   [cage C]

   entrance_antenna1 = 4

   entrance_antenna2 = 5

   
   [cage D]

   entrance_antenna1 = 6

   entrance_antenna2 = 7

   
   [cage A]

   entrance_antenna1 = 8

   entrance_antenna2 = 1


   [tunnel 1]

   entrance_antenna1 = 1

   entrance_antenna2 = 2


   [tunnel 2]

   entrance_antenna1 = 3

   entrance_antenna2 = 4


   [tunnel 3]

   entrance_antenna1 = 5

   entrance_antenna2 = 6


   [tunnel 4]

   entrance_antenna1 = 7

   entrance_antenna2 = 8

2. setup.txt file for the setup with internal antennas only. This file should be placed in the data directory with registrations of Eco-HAB setup with internal antennas.
   
   [setup]

   name = internal

   [cage A]

   internal_antenna = 1

   [cage C]

   internal_antenna = 8

3. Setup config file for the entire experiment:

   [shared compartment 1]

   setup_1_name = default
   
   compartment_1_name = cage A

   setup_2_name = internal
   
   compartment_2_name = cage A
   
   destination_name = cage A

   [shared compartment 2]

   setup_1_name = default
   
   compartment_1_name = cage C
   
   setup_2_name = internal
   
   compartment_2_name = cage C
   
   destination_name = cage C
 

   [rename compartment 1]

   setup_name = default
   
   compartment_name = cage B
   
   destination_name = cage B

   
   [rename compartment 2]
   
   setup_name = default
   
   compartment_name = cage D
   
   destination_name = cage D

This config file consists of two parts. The first part consisting of sections [shared compartment 1] and [shared compartment 2] specifies parts of the experimental setups that are shared by both submodules. In this case it is cage A, which has two entrance antennas, which are part of the setup named default, and an entrance antenna, which is a part of the setup named internal, and cage C. In this sections we specify locations and set the name that will be used in results files (in this case cage A and cage C). For clarity pyEcoHAB, when merging different setups into one modular dataset, adds setup names to names of the cages and tunnels that are not shared by different setups. One can rename these locations for easier further data analysis.


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



