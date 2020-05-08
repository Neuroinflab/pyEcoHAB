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
   
An Eco-HAB "experiment" is described by a :class:`pyEcoHAB.Loader`
object.  As an example, let us load data from
pyEcoHAB.sample_data. This path contains data from an experiment on
valproate-treated (VPA) BALB/c mice, which is a widely used strain.

Load the experiment:

>>> data = pyEcoHAB.Loader(pyEcoHAB.sample_data)

Eco-HAB data is stored as text files containing hourly registrations of animals by antennas. Each event is a separate line, in which id, date, time of registration, antenna id, duration of registration (in ms) and animal tag. :class:`pyEcoHAB.Loader` locates these hourly files in provided location (in this case it is pyEcoHAB.sample_data), reads them in and parses the events. :class:`pyEcoHAB.Loader` also performs data diagnostics:
1. Eliminates events that are not trustworthy ('ghost-tags').
2. Lists out which antennas had no registrations for a specified amount of time. Typically antennas that were silent for more that 1 h are listed together with date and time they were silent.
3. Counts how many antenna registrations have been omitted,


Timeline of the experiment 



