import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
     name='pyEcoHAB',
     version='1.0',
     author="Joanna Jędrzejewska-Szmek, Jan Mąka, Szymon Łęski",
     author_email="j.jedrzejewska-szmek@nencki.edu.pl",
     description="Read in and analyse mice behavioral data acquired by Eco-HAB",
     url="https://github.com/Neuroinflab/pyEcoHAB",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: BSD",
         "Operating System :: OS Independent",
     ],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    project_urls={
        "Bug Tracker": "https://github.com/Neuroinflab/pyEcoHAB/issues",
    },
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    include_package_data=True,
    package_data={'pyEcoHAB': ['data/*',
                               'data/BALB_VPA_data_cohort_1/*',
                               'data/BALB_VPA_data_cohort_1_divided/*',
                               'data/BALB_VPA_data_cohort_1_divided/*/*',
                               'data/empty',
                               'data/modular_1/*',
                               'data/modular_1/*/*',
                               'data/test_experiment_setups/*',
                               'data/test_setups/*',
                               'data/test_setups_2/*',
                               'data/time_change/*',
                               'data/weird_3_mice/*',
                               'data/weird_short/*',
                               'data/weird_short_3_mice/*',
                               'data/weird_very_short/*',
                               'data/weird_very_short_3_mice/*',
                               ]},
    license='LGPL-2.1-or-later',
 )
