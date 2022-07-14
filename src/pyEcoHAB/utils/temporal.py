from datetime import datetime, timedelta
from pyEcoHAB.utils import for_loading as fl


def convert_int_to_time(number):
    if number < 10:
        return "0%d" % number
    return "%d" % number


def find_light_beginning(dark_beg, dark_len):
    hours, mins = dark_beg.split(":")
    d_hours = int(hours)
    d_mins = int(mins)
    # dark len in h
    add_to_h = 0
    new_mins = int((dark_len - int(dark_len))*60) + d_mins
    if new_mins > 60:
        new_mins = new_mins-60
        add_to_h = 1
    new_hours = d_hours + int(dark_len) + add_to_h
    if new_hours >= 24:
        new_hours = new_hours - 24
    return "%s:%s" % (convert_int_to_time(new_hours),
                      convert_int_to_time(new_mins))


def find_first_last(filename_list):
    first = filename_list[0].split("_")[0]
    last = filename_list[-1].split(".txt")[0] + " UTC"
    return first, last


def last_day_to_datetime(last_day):
    return datetime.strptime(last_day, "%Y%m%d_%H%M%S %Z")


def strtime_to_datetime(text_time):
    return datetime.strptime(text_time, "%Y%m%d%H:%M %Z")


def get_date(date_obj):
    return date_obj.strftime("%d.%m.%Y")


def get_time(date_obj):
    return date_obj.strftime("%H:%M")


def make_config_entry(start_date, end_date):
    return {
        "startdate": get_date(start_date),
        "starttime": get_time(start_date),
        "enddate": get_date(end_date),
        "endtime": get_time(end_date),
        }


def gen_timeline(data_directory, dark_beginning="12:00",
                 first_phase="dark", dark_length=12,
                 light_length=12,
                 phase_name="EMPTY"):
    """
    Automatically generate timeline of an EcoHAB experiment

    This function will calculate phases beginnings and endings and generate
    the timeline of the experiment with phases named phase_name number
    phase_type. The file will be save in data_directory. If the beginning of
    the dark phase is not provided, 12:00 will be used. Dark and light phase
    lengths will be used to calculate begginings and endings of each phase.
    Phase lengths can be specified, otherwise it is assumed that dark and light
    phase are 12 h long.


    Args:
       data_directory: str
           path to the directory containing experiment data files
       dark_beginning: str
           At what time do the dark phases begin. Default: 12:00
       first_phase: str
           What phase is the first: dark or light. Default: dark
       dark_length: float
           Length of the dark phase (in hours). Default: 12
       light_length: float
           Length of the light phase (in hours). Default: 12
       phase_name: str
           name of all the phases. Default: EMPTY.
           Consecutive phases will be named: EMPTY 1 dark, EMPTY 1 light,
           EMPTY 2 dark ...
    """
    config = {}
    # find files
    filenames = sorted(fl.get_filenames(data_directory))
    # find beginning of the experiment
    first_day, last_day = find_first_last(filenames)
    light_beginning = find_light_beginning(dark_beginning,
                                           dark_length)
    light_duration = timedelta(hours=light_length)
    dark_duration = timedelta(hours=dark_length)

    if first_phase.lower() == "dark":
        str_date = "%s%s UTC" % (first_day, dark_beginning)

    elif first_phase.lower() == "light":
        str_date = "%s%s UTC" % (first_day, light_beginning)

    start_date = strtime_to_datetime(str_date)
    total_beg = strtime_to_datetime(str_date)
    total_end = last_day_to_datetime(last_day)
    i = 1
    current_phase = first_phase
    while True:
        # current phase name
        full_phase_name = "%s %d %s" % (phase_name, i, current_phase)
        if current_phase.lower() == "light":
            end_date = start_date + light_duration
        else:
            end_date = start_date + dark_duration
        config[full_phase_name] = make_config_entry(start_date,
                                                    end_date)
        if current_phase.lower() == "light":
            i = i+1
            current_phase = "dark"
        else:
            current_phase = "light"
        start_date = end_date
        if start_date > total_end:
            config["ALL"] = make_config_entry(total_beg,
                                              end_date)
            break
    return config
