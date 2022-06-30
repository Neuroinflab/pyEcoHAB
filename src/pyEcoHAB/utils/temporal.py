import time

def convert_int_to_time(number):
    if number < 10:
        return "0%d"%number
    return "%d"%number

def find_light_beginning(dark_beg, dark_len):
    hours, mins = dark_beg.split(":")
    d_hours = int(hours)
    d_mins = int(mins)
    #dark len in h
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
    return filename_list[0].split("_")[0], filenames[-1].split("_"][0]
