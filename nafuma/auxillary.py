import json
import numpy as np
import os

def update_options(options, required_options, default_options):
    ''' Takes a dictionary of options along with a list of required options and dictionary of default options, and sets all keyval-pairs of options that is not already defined to the default values'''

    for option in required_options:
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
    from datetime import datetime

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