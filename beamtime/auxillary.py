import json

def update_options(options, required_options, default_options):
    ''' Takes a dictionary of options along with a list of required options and dictionary of default options, and sets all keyval-pairs of options that is not already defined to the default values'''

    for option in required_options:
        if option not in options.keys():
            options[option] = default_options[option]


    return options

def save_options(options, path):
    ''' Saves any options dictionary to a JSON-file in the specified path'''

    with open(path, 'w') as f:
        json.dump(options,f)


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



def hello_world2(a=1, b=2):

    print(f'Halla, MAFAKKAS! a = {a} og b = {b}')