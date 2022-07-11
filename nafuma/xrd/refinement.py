import os
import shutil
import time
import datetime

import nafuma.auxillary as aux

def make_big_inp(data: dict, options={}):

    required_options = ['output', 'overwrite', 'backup', 'backup_path']

    default_options = {
        'output': 'big.inp', # Name of the output .INP file
        'overwrite': False, # Toggles overwrite on / off
        'backup': True, # Toggles backup on / off. Makes a backup of the file if it already exists. Only runs if overwrite is enabled.
        'backup_dir': 'backup' # Specifies the path where the backup files should be located
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    # Raises exception if files exists and overwrite is not enabled
    if not options['overwrite']:
        if os.path.exists(options['output']):
            raise Exception(f'Overwrite disabled and file already exists: {options["output"]}')


    # Makes a backup of file
    elif options['backup'] and os.path.exists(options['output']):
        aux.backup_file(filename=options['output'], backup_dir=options['backup_dir'])



    