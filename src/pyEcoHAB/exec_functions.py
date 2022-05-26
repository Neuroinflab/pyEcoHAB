# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
from .write_to_file import save_single_histograms, write_csv_rasters
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .plotting_functions import single_heat_map
from .utils import general as utils


def evaluate_whole_experiment(ecohab_data, timeline, main_directory,
                              prefix, func, fname,
                              xlabel, ylabel, title,
                              args=[], remove_mouse=None,
                              vmin=None, vmax=None,
                              delimiter=";"):

    phases = utils.filter_dark_light(timeline.sections())
    mice = [mouse[-4:] for mouse in ecohab_data.mice]
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    result = np.zeros((len(phases), len(mice), len(mice)))
    fname_ = '%s_%s%s.csv' % (fname, prefix, add_info_mice)
    hist_dir = os.path.join("other_variables", fname, 'histograms')
    rast_dir = os.path.join("other_variables", fname, 'raster_plots')
    for i, phase in enumerate(phases):
        if len(args):
            result[i] = func(ecohab_data, timeline, phase, *args)
        else:
            result[i] = func(ecohab_data, timeline, phase)
        save_single_histograms(result[i],
                               fname,
                               ecohab_data.mice,
                               phase,
                               main_directory,
                               hist_dir,
                               prefix,
                               additional_info=add_info_mice,
                               delimiter=delimiter)
        single_heat_map(result[i],
                        fname,
                        main_directory,
                        mice,
                        prefix,
                        phase,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        subdirectory=hist_dir,
                        vmax=vmin,
                        vmin=vmax,
                        xticks=mice,
                        yticks=mice)
    write_csv_rasters(ecohab_data.mice,
                      phases,
                      result,
                      main_directory,
                      rast_dir,
                      fname_,
                      delimiter=delimiter,
                      symmetrical=False, prefix=prefix)
    make_RasterPlot(main_directory,
                    rast_dir,
                    result,
                    phases,
                    fname_,
                    ecohab_data.mice,
                    title=title,
                    symmetrical=False, prefix=prefix)
