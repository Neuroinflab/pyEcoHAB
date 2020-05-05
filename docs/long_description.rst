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

Eco-HAB data is stored in text files containing hourly registrations of animal tags by nearby antennas. In the files each event is a separate line containing registration id, date, time of registration, antenna id, duration of registration (in ms) and animal tag. :class:`pyEcoHAB.Loader` locates these hourly files in provided location (in this case it is pyEcoHAB.sample_data), reads them in and parses the events. :class:`pyEcoHAB.Loader` also performs data diagnostics:
1. Eliminates events that are not trustworthy ('ghost-tags').
2. Lists out which antennas had no registrations for a specified amount of time. Typically antennas that were silent for more that 1 h are listed together with date and time they were silent.
3. Counts how many antenna registrations have been omitted. The design of the Eco-HAB restricts, which antennas recorded an animal tag in a sequence. For a specified animal tag for consecutive recordings antenna id can differ either by 1 or by 7 (in case of antennas 1 and 8). :class:`pyEcoHAB.Loader` counts all the cases when antenna id difference is incorrect and lists the statistics as in-correct registrations.


Timeline of the experiment 



