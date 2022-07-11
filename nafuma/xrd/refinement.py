import os
import shutil
import time
import re
import datetime
import warnings

import nafuma.auxillary as aux

def make_big_inp(data: dict, options={}):
    ''' Generates a big .INP-file with all filenames found in data["path"]. Uses a template .INP-file (which has to be generated manually from an initial refinement in TOPAS) and appends this to a large .INP-file
    while changing the filenames. '''

    required_options = ['template', 'output', 'overwrite', 'backup', 'backup_dir', 'include', 'topas_options', 'save_dir', 'log', 'logfile']

    default_options = {
        'template': 'start.inp', # Name of the template .INP-file
        'output': 'big.inp', # Name of the output .INP-file
        'overwrite': False, # Toggles overwrite on / off
        'backup': True, # Toggles backup on / off. Makes a backup of the file if it already exists. Only runs if overwrite is enabled.
        'backup_dir': 'backup', # Specifies the path where the backup files should be located
        'include': [],
        'topas_options': {
            'bootstrap_errors': None,
            'A_matrix_memory_allowed_in_Mbytes': None,
            'approximate_A': False,
            'conserve_memory': False,
            'do_errors': False,
            'continue_after_convergence': False,
            'num_runs': None,
        },
        'save_dir': 'results',
        'log': False,
        'logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_generate_big_inp.log',
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if not os.path.exists(options['template']):
        raise Exception(f'Template file not found: {options["template"]}')

    # Raises exception if files exists and overwrite is not enabled
    if not options['overwrite']:
        if os.path.exists(options['output']):
            raise Exception(f'Overwrite disabled and file already exists: {options["output"]}')


    # Makes a backup of file
    elif options['backup'] and os.path.exists(options['output']):

        if options['log']:
            aux.write_log(message=f'File {options["output"]} already exists. Creating backup in {options["backup_dir"]}.')

        aux.backup_file(filename=options['output'], backup_dir=options['backup_dir'])



    runlist = os.path.join(os.path.dirname(options['output']), 'runlist.txt')
    options['include'].append(runlist)

    with open(options['template'], 'r') as template, open(options['output'], 'w', newline='\n') as output, open(runlist, 'w', newline='\n') as runlist:
             
        write_headers(output, options)
        template = template.read()
        
        for i, path in enumerate(data['path']):

            s = change_labels_and_paths(template, path, i, options)
            output.write(s)
            
            runlist.write('#define \tUSE_'+f'{i}'.zfill(4) + '\n')



def write_headers(fout, options):

    for file in options['include']:
        fout.write(f'#include {file} \n')

    fout.write('\n')

    for option, value in options['topas_options'].items():
        if value:
            fout.write(f'{option} {value} \n')



def change_labels_and_paths(template, path, i, options):

    temp_xdd = find_xdd_in_inp(options['template'])
    
    # Replace diffractogram-path
    s = template.replace(temp_xdd, path).replace('XXXX', f'{i}'.zfill(4))

    basename = os.path.basename(path).split(".")[0]
    
    # Define regular expressions for output lines
    regs = [r'Out_Riet\([\S]*\)',
            r'Out_CIF_STR\([\S]*\)',
            r'Out_CIF_ADPs\([\S]*\)',
            r'Out_CIF_Bonds_Angles\([\S]*\)',
            r'Out_FCF\([\S]*\)',
            r'Create_hklm_d_Th2_Ip_file\([\S]*\)',
            r'out(.*?)append']
    
    # Define substitute strings for output lines
    subs = [f'Out_Riet({options["save_dir"]}/{basename}_riet.xy)',
            f'Out_CIF_STR({options["save_dir"]}/{basename}.cif))',
            f'Out_CIF_ADPs({options["save_dir"]}/{basename}.cif))',
            f'Out_CIF_Bonds_Angles({options["save_dir"]}/{basename}.cif))',
            f'Out_FCF({options["save_dir"]}/{basename}.fcf)',
            f'Create_hklm_d_Th2_Ip_file({options["save_dir"]}/{basename}_hkl.dat)',
            f'out \t {options["save_dir"]}/{basename}_refined_params.dat']

    # Substitute strings in output lines
    for reg, sub in zip(regs, subs):
        s = re.sub(reg, sub, s)


    return s



def find_xdd_in_inp(path):
    ''' Finds the path to the .xy / .xye scan in a given .INP-file. Assumes only one occurence of xdd and will return the last one no matter what, but outputs a UserWarning if more than one xdd is found.'''

    with open(path, 'r') as f:
        lines = f.readlines()


    xdds = 0
    for line in lines:
        if 'xdd' in line:
            xdd = line.split()[-1].strip('"')
            xdds += 1

    if xdds > 1:
        warnings.warn(f'More than one path was found in {path}. Returning last occurence - please make sure this is what you want!')

    return xdd
