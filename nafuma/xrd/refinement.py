from email.policy import default
import os
import shutil
import subprocess
import re

import time
import datetime

import warnings

import nafuma.auxillary as aux


def make_big_inp(data: dict, options={}):
    ''' Generates a big .INP-file with all filenames found in data["path"]. Uses a template .INP-file (which has to be generated manually from an initial refinement in TOPAS) and appends this to a large .INP-file
    while changing the filenames. '''

    required_options = ['template', 'output', 'overwrite', 'backup', 'backup_dir', 'include', 'topas_options', 'save_results', 'save_dir', 'log', 'logfile']

    default_options = {
        'template': 'start.inp', # Name of the template .INP-file
        'output': 'big.inp', # Name of the output .INP-file
        'overwrite': False, # Toggles overwrite on / off
        'backup': True, # Toggles backup on / off. Makes a backup of the file if it already exists. Only runs if overwrite is enabled.
        'backup_dir': 'backup', # Specifies the path where the backup files should be located
        'include': [],
        'topas_options': None,
        'save_results': True,
        'save_dir': 'results',
        'log': False,
        'logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_generate_big_inp.log',
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    required_topas_options = ['iters', 'chi2_convergence_criteria', 'bootstrap_errors', 'A_matrix_memory_allowed_in_Mbytes', 'approximate_A', 'conserve_memory', 'do_errors', 'continue_after_convergence', 'num_runs']
    default_topas_options = {
            'iters': 100000,
            'chi2_convergence_criteria': 0.001,
            'bootstrap_errors': None,
            'A_matrix_memory_allowed_in_Mbytes': None,
            'approximate_A': False,
            'conserve_memory': False,
            'do_errors': False,
            'continue_after_convergence': False,
            'num_runs': None,
        }

    options['topas_options'] = aux.update_options(options=options['topas_options'], required_options=required_topas_options, default_options=default_topas_options)


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

            s = make_inp_entry(template=template, xdd=path, num=i, options=options)
            output.write(s)
            output.write('\n\n')
            
            runlist.write('#define \tUSE_'+f'{i}'.zfill(4) + '\n')



def write_headers(fout, options):
    # FIXME Could modify this to make sure that certain combinations of options is not written, such as both do_errors and bootstrap_errors. 

    for file in options['include']:
        fout.write(f'#include {file} \n')

    fout.write('\n')

    if options['save_results']:
        fout.write('#define output \n')
    
    fout.write('\n')

    for option, value in options['topas_options'].items():
        if value:
            if isinstance(value, bool):
                fout.write(f'{option} \n')
            else:
                fout.write(f'{option} {value} \n')


def get_headers(inp):

    with open(inp, 'r') as inp:
        headers = []

        line = inp.readline()

        while not all(keyword in line for keyword in ['out', 'append']):
            line = inp.readline()

        # Jump down to lines
        line = inp.readline()
        line = inp.readline()

        while not 'Out_String' in line:
                
            if line.split():
                header = line.split()[1]
                if all(keyword in header for keyword in ['Get', '(', ')']):
                    header = header[4:-1]

                headers.append(header)

            line = inp.readline()


    return headers
            

def get_paths(inp):

    paths = []

    with open(inp, 'r') as inp:
        lines = inp.readlines()

    for line in lines:
        if all(keyword in line for keyword in ['out', 'append']):
            paths.append(line.split()[1])


    return paths

def make_inp_entry(template: str, xdd: str, num: int, options: dict) -> str:
    ''' Takes a template and creates an entry with xdd as path to file and with number num.'''

    temp_xdd = find_xdd_in_inp(options['template'])


    num_str = f'{num}'.zfill(4)
    
    # Replace diffractogram-path
    s = template.replace(temp_xdd, xdd).replace('XXXX', num_str)

    basename = os.path.basename(xdd).split(".")[0]
    
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
            f'Out_CIF_STR({options["save_dir"]}/{basename}.cif)',
            f'Out_CIF_ADPs({options["save_dir"]}/{basename}.cif)',
            f'Out_CIF_Bonds_Angles({options["save_dir"]}/{basename}.cif)',
            f'Out_FCF({options["save_dir"]}/{basename}.fcf)',
            f'Create_hklm_d_Th2_Ip_file({options["save_dir"]}/{basename}_hkl.dat)',
            f'out \t {options["save_dir"]}/{basename}_refined_params.dat \t append']

    # Substitute strings in output lines
    for reg, sub in zip(regs, subs):
        s = re.sub(reg, sub, s)



    return s



def find_xdd_in_inp(path: str) -> str:
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


def refine(data: dict, options={}):
    ''' Calls TOPAS from options['topas_path'], which should point to tc.exe in the TOPAS-folder. If not explicitly passed, will try to use what is in topas.conf.'''

    required_options = ['topas_path', 'topas_log', 'topas_logfile', 'log', 'logfile', 'overwrite']

    default_options = {
        'topas_path': None,
        'topas_log': False,
        'topas_logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_topas.out',
        'log': False,
        'logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_generate_big_inp.log',
        'overwrite': False
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if not options['topas_path']:
        import nafuma.xrd as xrd

        # Open topas.conf in the nafuma/xrd 
        with open(os.path.join(os.path.dirname(xrd.__file__), 'topas.conf'), 'r') as f:
            topas_base = f.read()

        options['topas_path'] = os.path.join(topas_base, 'tc.exe')


    # Check to see if the executable exists
    if not os.path.exists(options['topas_path']):
        raise Exception('TOPAS executable not found! Please explicitly pass path to tc.exe directly in options["topas_path"] or change base folder in topas.conf to the correct one.')


    # Create folders if they don't exist
    paths, headers = get_paths(data['inp']), get_headers(data['inp'])
    
    for path in paths:
        dirname = os.path.dirname(path)

        if dirname and not os.path.isdir(dirname):
            os.makedirs(dirname)

        if not os.path.exists(path) or options['overwrite']:
            with open(path, 'w') as results:
                for header in headers:
                    results.write(header+'\t')

                results.write('\n')
        else:
            raise Exception(f'Results file already exists: {path}')


    # Create shell command
    command = ' '.join([options['topas_path'], data['inp']])

    # Append output if logging is enabled
    if options['topas_log']:
        command = ' '.join([command, f'>{options["topas_logfile"]}'])

        if os.path.dirname(options['topas_logfile']) and not os.path.isdir(os.path.dirname(options['topas_logfile'])):
            os.makedirs(os.path.dirname(options['topas_logfile']))

    
    subprocess.call(command, shell=True)