How to use :doc:`pyEcoHAB` to load and analyze an :doc:`Eco-HAB` datatest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Loading and analyzing data consists of the following steps:

.. contents::
   :local:

The short version is:

#install numpy and matplotlib

#clone pyEcoHAB

#get into python using:

>PYTHONPATH=$PYTHONPATH:/full/path/to/pyEcoHAB/: python3

>>> import pyEcoHAB
>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)
>>> config = pyEcoHAB.ExperimentConfigFile(pyEcoHAB.sample_data)
>>> pyEcoHAB.get_activity(data, config, 3600)
>>> pyEcoHAB.get_incohort_sociability(data, config, 3600)
>>> pyEcoHAB.get_solitude(data, config)
>>> pyEcoHAB.get_following(data, config, 1000)

The long version is below.

Experimental data
``````````````````````
   
An Eco-HAB "experiment" is described by a :class:`pyEcoHAB.Loader` object.
As an example, let's load data from pyEcoHAB.sample_data. This path contains data
from an experiment on valproate-treated (VPA) BALB/c mice, which is a widely used strain.

Load the experiment:

>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)

Timeline of the experiment 



