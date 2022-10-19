import json
import numpy as np
import os
import shutil

import time
from datetime import datetime

def update_options(options, default_options, required_options=None):
    ''' Takes a dictionary of options along with a list of required options and dictionary of default options, and sets all keyval-pairs of options that is not already defined to the default values'''
    #FIXME This has been updated so that required_options is not needed. But lots of scripts still passes required_options, so for now it is still accepted, but has a default value and remains unused. Needs to go through all scripts to stop passing of this variable to remove it.


    for option in default_options.keys():
        if option not in options.keys():
            options[option] = default_options[option]


    return options

def save_options(options, path, ignore=None):
    ''' Saves any options dictionary to a JSON-file in the specified path'''

    options_copy = options.copy()

    if ignore:
        if not isinstance(ignore, list):
            ignore = [ignore]

        for i in ignore:
            options_copy[i] = 'Removed'


    if not os.path.isdir(os.path.dirname(path)):
        if os.path.dirname(path):
            os.makedirs(os.path.dirname(path))


    with open(path, 'w') as f:
        json.dump(options_copy, f, skipkeys=True, indent=4)


def load_options(path):
    ''' Loads JSON-file into a dictionary'''

    with open(path, 'r') as f:
        options = json.load(f)

    return(options)



def swap_values(dict, key1, key2):

    key1_val = dict[key1]
    dict[key1] = dict[key2]
    dict[key2] = key1_val

    return dict



def ceil(a, roundto=1):

    fac = 1/roundto

    a = np.ceil(a*fac) / fac

    return a

def floor(a, roundto=1):

    fac = 1/roundto

    a = np.floor(a*fac) / fac

    return a



def write_log(message, options={}):


    required_options = ['logfile']
    default_options = {
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S.log")}'
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    if not os.path.isdir(os.path.dirname(options['logfile'])):
        os.makedirs(os.path.dirname(options['logfile']))


    now = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    message = f'[{now}] {message} \n'


    with open(options['logfile'], 'a') as f:
        f.write(message)


#Function that "collects" all the files in a folder, only accepting .dat-files from xanes-measurements
def get_filenames(path, ext, filter=''):
    ''' Collects all filenames from specified path with a specificed extension
    
    Input:
    path: path to find all filenames (relative or absolute)
    ext: extension (including ".")'''
 
    filenames = [os.path.join(path, filename) for filename in os.listdir(path) if os.path.isfile(os.path.join(path, filename)) and filename.endswith(ext) and filter in filename] 
        
    return filenames

def move_list_element_last(filenames,string):
    for i,file in enumerate(filenames):
        if string in file:
            del filenames[i]
            filenames.append(file)
    return filenames



def backup_file(filename, backup_dir):
    # Creates backup-folder if it does not exist
    if not os.path.isdir(backup_dir):
        os.makedirs(backup_dir)
    

    # Get a list of all previous backup files with the same basename as well as the creation time for the 
    prev_backup_files = [file for file in os.listdir(backup_dir) if os.path.basename(filename.split('.')[0]) in file]
    creation_time = datetime.strptime(time.ctime(os.path.getmtime(filename)), '%a %b %d %H:%M:%S %Y').strftime("%Y-%m-%d_%H-%M-%S")
    ext = '.' + filename.split('.')[-1]

    dst_basename = creation_time + '_' + filename.split('.')[0] + '_' + f'{len(prev_backup_files)}'.zfill(4) + ext
    dst = os.path.join(backup_dir, dst_basename)


    shutil.copy(filename, dst)


def get_unique(full_list):

    unique_list = []

    for entry in full_list:
        if not entry in unique_list:
            unique_list.append(entry)

    return unique_list


def swap_values(options: dict, key1, key2):

    if not isinstance(key1,list):
        key1 = [key1]
    if not isinstance(key2,list):
        key2 = [key2]

    assert len(key1) == len(key2)

    for k1, k2 in zip(key1, key2):
        options[k1], options[k2] = options[k2], options[k1]


    return options



def find_neighbours(value, df, colname, start=0, end=-1):
    ''' Finds closest match to a given value in colname of df. If there is an exact match, returns index of this value. Else, it returns the nearest neighbors (upper and lower)'''
    
    df = df.iloc[start:end]
    
    exactmatch = df[df[colname] == value]
    if not exactmatch.empty:
        return exactmatch.index
    else:
        lower_df = df[df[colname] < value][colname]
        upper_df = df[df[colname] > value][colname]

        lowerneighbour_ind = lower_df.idxmax()
        upperneighbour_ind = upper_df.idxmin()

        return [lowerneighbour_ind, upperneighbour_ind] 