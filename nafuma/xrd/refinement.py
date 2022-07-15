import os
import shutil
import subprocess
import re

import time
import datetime
import warnings
import json

import pandas as pd

import nafuma.auxillary as aux



def make_initial_inp(data: dict, options={}):

    required_options = ['filename', 'overwrite', 'include', 'save_results', 'save_dir', 'topas_options', 'background', 'capillary', 'start', 'finish', 'exclude']

    default_options = {
        'filename': 'start.inp',
        'overwrite': False,
        'save_results': False,
        'save_dir': 'results/',
        'include': [], # Any files that should be included with #include
        'topas_options': {}, 
        'background': 7,
        'capillary': False,
        'interval': [None, None], # Start and finish values that TOPAS should refine on. Overrides 'start' and 'finish'
        'start': None, # Start value only. Overridden by 'interval' if this is set.
        'finish': None, # Finish value only. Overridden by 'interval' if this is set.
        'exlude': [] # Excluded regions. List of lists.
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    options['topas_options'] = update_topas_options(options=options)


    if not os.path.exists(options['filename']) or options['overwrite']:
        with open(options['filename'], 'w') as fout:
            write_headers(fout=fout, options=options)


            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
            fout.write('\'START OF INP - XXXX')
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
            fout.write('r_wp    0    r_exp    0    r_p    0    r_wp_dash    0    r_p_dash    0    r_exp_dash    0    weighted_Durbin_Watson    0    gof    0')
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')

            write_xdd(fout=fout, data=data, options=options)

            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
            fout.write('\'PARAMETER DEFINITONS')
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')

            for i, str in enumerate(data['str']):
                fout.write('\n\n')
                write_params(fout=fout, data=data, options=options, index=i)


            
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
            fout.write('\'STR DEFINITIONS')
            fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')

            for i, str in enumerate(data['str']):
                fout.write('\n\n')
                write_str(fout=fout, data=data, options=options, index=i)


           

def write_xdd(fout, data, options):
    import nafuma.xrd as xrd
    basedir = os.path.dirname(xrd.__file__)

    with open (os.path.join(basedir, 'snippets.json'), 'r') as snip:
        snippets = json.load(snip)

    topas_options = options['topas_options']
    
    # Write initial parameters
    fout.write(f'xdd "{data["xdd"]}"\n')
    fout.write('\t'+snippets['calculation_step'].format(topas_options['convolution_step'])+'\n')
    
    for include in options['include']:
        if 'peak_width' in include:
            with open(include, 'r') as f:
                lines = f.readlines()
                peak_width_params = []
                for line in lines:
                    param = line.split()[1][1:]
                    peak_width_params.append(param)



            fout.write('\t'+snippets['gauss_fwhm'].format(peak_width_params[0], peak_width_params[1], peak_width_params[2])+'\n')
    
    # Write background
    fout.write('\tbkg @    ')
    for i in range(options['background']):
        fout.write('0    ')

    # Write wavelength and LP-factor
    fout.write('\n\n')
    
    fout.write('\t'+snippets['wavelength'].format(data['wavelength'])+'\n')
    fout.write('\t'+snippets['lp_factor'].format(topas_options['lp_factor'])+'\n')
    fout.write('\t'+snippets['zero_error']+'\n')

    if options['capillary']:
        fout.write('\n')
        for i, line in enumerate(snippets['capillary']):
            if i == 0:
                line = line.format(topas_options['packing_density'])
            if i == 1:
                line = line.format(topas_options['capdia'])

            fout.write('\t'+line+'\n')


    fout.write('\n')
    if options['interval'][0] or options['start']:
        options['start'] = options['interval'][0] if options['interval'][0] else options['start']
        fout.write(f'\tstart_X {options["start"]}\n')

    if options['interval'][1] or options['finish']:
        options['finish'] = options['interval'][1] if options['interval'][1] else options['finish']
        fout.write(f'\tfinish_X {options["finish"]}\n')

    if options['exclude']:
        for exclude in options['exclude']:
            fout.write(f'\texclude {exclude[0]} {exclude[1]}\n')




def write_params(fout, data, options, index=0):

    atoms = read_cif(data['str'][index])
    if 'labels' in data.keys():
        label = data['labels'][index]
    else:
        label = index


    fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
    fout.write(atoms['_chemical_name_common'] + f'({atoms["_space_group_name_H-M_alt"]}) - Parameters')
    fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')

    a       = float(atoms['_cell_length_a'].split('(')[0])
    b       = float(atoms['_cell_length_b'].split('(')[0])
    c       = float(atoms['_cell_length_c'].split('(')[0])
    alpha   = float(atoms['_cell_angle_alpha'].split('(')[0])
    beta    = float(atoms['_cell_angle_beta'].split('(')[0])
    gamma   = float(atoms['_cell_angle_gamma'].split('(')[0])


    # WRITE LATTICE PARAMETERS
    # If start_values is defined:
    fout.write(f'#ifdef start_values_{label}\n')
    lpa = f'local !lpa_{label} {a} ;: {a}'
    lpb = f'local !lpb_{label} {b} ;: {b}'
    lpc = f'local !lpc_{label} {c} ;: {c}'
    fout.write('{: <55}  {: <55}  {: <55}\n'.format(lpa, lpb, lpc))
    lpal = f'local !lpal_{label} {alpha} ;: {alpha}'
    lpbe = f'local !lpbe_{label} {beta} ;: {beta}'
    lpga = f'local !lpga_{label} {gamma} ;: {gamma}'
    fout.write('{: <55}  {: <55}  {: <55}\n'.format(lpal, lpbe, lpga))

    # Otherwise
    fout.write('\n')
    fout.write('#else\n')
    lpa = f'local !lpa_{label} {a}'
    lpb = f'local !lpb_{label} {b}'
    lpc = f'local !lpc_{label} {c}'
    fout.write('{: <55}  {: <55}  {: <55}\n'.format(lpa, lpb, lpc))
    lpal = f'local !lpal_{label} {alpha}'
    lpbe = f'local !lpbe_{label} {beta}'
    lpga = f'local !lpga_{label} {gamma}'
    fout.write('{: <55}  {: <55}  {: <55}\n'.format(lpal, lpbe, lpga))
    fout.write('#endif\n\n')

    sites = list(atoms['atoms'].keys())

    attrs = {
        '_atom_site_fract_x': 'x', 
        '_atom_site_fract_y': 'y', 
        '_atom_site_fract_z': 'z', 
        '_atom_site_occupancy': 'occ',
        '_atom_site_B_iso_or_equiv': 'beq'
        }


    # WRITE SITE PARAMETERS
    for site in sites:

        params = []
        for attr in attrs:
            if attr in atoms["atoms"][site].keys():
                value = atoms["atoms"][site][attr].split("(")[0]
                value = value if value != '.' else 0.

                params.append('{: <20} {: <20}'.format(f'local !{attrs[attr]}_{site}_{label}', f' = {value} ;: {value}'))
                #fout.write(f'local {attrs[attr]}_{site}_{label}\t\t =\t {value} \t ;= \t\t\t\t\t')

        fout.write('{: <55}  {: <55}  {: <55} {: <55} {: <55}\n'.format(*params))

    fout.write('\n')

    fout.write('{: <55}  {: <55}  {: <55} {: <55}\n'.format(
                                                        f'local !csgc_{label}_XXXX = 200 ;: 200',
                                                        f'local !cslc_{label}_XXXX = 200 ;: 200',
                                                        f'local !sgc_{label}_XXXX = 0 ;: 0',
                                                        f'local !slc_{label}_XXXX = 0 ;: 0',
    ))


    # fout.write(f'local csgc_{label} ;:\t\t\t\t\t')
    # fout.write(f'local cslc_{label} ;:\t\t\t\t\t')
    # fout.write(f'local sgc_{label} ;:\t\t\t\t\t')
    # fout.write(f'local slc_{label} ;:\n')




def write_str(fout, data, options, index=0):

    atoms = read_cif(data['str'][index])
    if 'labels' in data.keys():
        label = data['labels'][index]
    else:
        label = index

   
    fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')
    fout.write('\'' + atoms['_chemical_name_common'].strip('\'') + ' ' + f'({atoms["_space_group_name_H-M_alt"]})'.replace("\'", ""))
    fout.write('\n\'------------------------------------------------------------------------------------------------------------------------------------------\n')


    fout.write('\tstr\n')
    fout.write(f'\t\tphase_name "{atoms["_chemical_name_common"]} ({atoms["_space_group_name_H-M_alt"]})"\n'.replace('\'', ''))
    fout.write(f'\t\tspace_group {atoms["_space_group_IT_number"]}\n')
    fout.write(f'\t\ta = lpa_{label} ;\n')
    fout.write(f'\t\tb = lpb_{label} ;\n')
    fout.write(f'\t\tc = lpc_{label} ;\n')
    fout.write(f'\t\tal = lpal_{label} ;\n')
    fout.write(f'\t\tbe = lpbe_{label} ;\n')
    fout.write(f'\t\tga = lpga_{label} ;\n')
    fout.write('\n')
    #FIXME fix the if-statement below, so that cell-volume is the correct formula, based on lattice params and angles.
    if '_cell_volume' not in atoms.keys():
        atoms['_cell_volume'] = 0
    else:
        fout.write(f'\t\tcell_volume\t vol_{label}_XXXX {atoms["_cell_volume"]}\n')
    fout.write(f'\t\tcell_mass\t mass_{label}_XXXX 1\n')
    fout.write(f'\t\tweight_percent\t wp_{label}_XXXX 100\n\n')
    fout.write('\n')
    fout.write('\t\tscale @ 1.0\n')
    fout.write('\t\tr_bragg 1.0\n\n')

    fout.write(f'\t\t#ifdef crystallite_size_gaussian_{label}\n')
    fout.write(f'\t\tCS_G(csgc_{label}_XXXX)\n')
    fout.write('\t\t#endif\n')
    fout.write(f'\t\t#ifdef crystallite_size_lorentzian_{label}\n')
    fout.write(f'\t\tCS_L(cslc_{label}_XXXX)\n')
    fout.write('\t\t#endif\n\n')

    fout.write(f'\t\t#ifdef strain_gaussian_{label}\n')
    fout.write(f'\t\tStrain_G(sgc_{label}_XXXX)\n')
    fout.write('\t\t#endif\n')
    fout.write(f'\t\t#ifdef strain_lorentzian_{label}\n')
    fout.write(f'\t\tStrain_L(slc_{label}_XXXX)\n')
    fout.write('\t\t#endif\n\n')





    for atom in atoms['atoms'].keys():

        specie = '{}{}'.format(
                                    atoms["atoms"][atom]["_atom_site_type_symbol"], 
                                    data['oxidation_states'][index][atoms["atoms"][atom]["_atom_site_type_symbol"]]
        )

        site = f'site {atom}'
        x = f'\t\tx = x_{atom}_{label}'
        y = f'\t\ty = y_{atom}_{label}'
        z = f'\t\tz = z_{atom}_{label}'
        occ = '\t\tocc {: <4} = occ_{}_{}'.format(specie, atom, label)
        #occ = f'\t\tocc {atoms["atoms"][atom]["_atom_site_type_symbol"]} occ_{atom}_{label}'
        beq = f'\t\tbeq = beq_{atom}_{label}'

        # FIXME Fix alignment here at some point
        #fout.write(f'\t\tsite {atom}\t x = x_{atom}_{label} ;\t y = y_{atom}_{label} ;\t z = z_{atom}_{label} ;\t occ {atoms["atoms"][atom]["_atom_site_type_symbol"]} occ_{atom}_{label} \t beq = beq_{atom}_{label} \t ;\n')
        fout.write('\t\t{: <9} {: <30}; {: <30}; {: <30}; {: <30}; {: <30};'.format(site, x, y, z, occ, beq))
        fout.write('\n')



    if options['save_results']:
        fout.write('\n\n')
        write_output(fout=fout, data=data, options=options, index=index)


    
def write_output(fout, data, options, index=0):


    filename = os.path.basename(data['xdd']).split('.')[0]

    atoms = read_cif(data['str'][index])

    if 'labels' in data.keys():
        label = data['labels'][index]
    else:
        label = index


    fout.write('#ifdef output\n')
    fout.write(f'\t\tOut_Riet({options["save_dir"]}/{label}_XXXX_riet.xy)\n')
    fout.write(f'\t\tOut_CIF_STR({options["save_dir"]}/{label}_XXXX.cif)\n')
    fout.write(f'\t\tOut_CIF_ADPs({options["save_dir"]}/{label}_XXXX.cif)\n')
    fout.write(f'\t\tOut_CIF_Bonds_Angles({options["save_dir"]}/{label}_XXXX.cif)\n')
    fout.write(f'\t\tCreate_hklm_d_Th2_Ip_file({options["save_dir"]}/{label}_XXXX_hkl.dat)\n')

    fout.write('\n')

    fout.write(f'out {options["save_dir"]}/{filename}_{label}.dat append\n')
    fout.write(f'\t\tOut_String("XXXX")\n')


    # FIXME Does not write out weighted_Durbin_Watson, TOPAS complained about this
    fout.write('\t\t{: <40} {: <40} {: <40} {: <40} {: <40} {: <40}'.format(
                                                                    f'Out(Get(r_wp), "%11.5f")',
                                                                    f'Out(Get(r_exp), "%11.5f")',
                                                                    f'Out(Get(r_p), "%11.5f")',
                                                                    f'Out(Get(r_p_dash), "%11.5f")',
                                                                    f'Out(Get(r_exp_dash), "%11.5f")',
                                                                    f'Out(Get(gof), "%11.5f")',
        )
    )

    fout.write('\n')


    fout.write('\t\t{: <40} {: <40} {: <40}'.format(
                                                                    f'Out(vol_{label}_XXXX, "%11.5f")',
                                                                    f'Out(mass_{label}_XXXX, "%11.5f")',
                                                                    f'Out(wp_{label}_XXXX, "%11.5f")',
    
        )
    )

    
    fout.write('\n')


    fout.write('\t\t{: <40} {: <40} {: <40} {: <40} {: <40} {: <40}'.format(
                                                                    f'Out(lpa_{label}, "%11.5f")',
                                                                    f'Out(lpb_{label}, "%11.5f")',
                                                                    f'Out(lpc_{label}, "%11.5f")',
                                                                    f'Out(lpal_{label}, "%11.5f")',
                                                                    f'Out(lpbe_{label}, "%11.5f")',
                                                                    f'Out(lpga_{label}, "%11.5f")',
    
        )
    )

    fout.write('\n\n')

    for atom in atoms['atoms']:
        fout.write('\t\t{: <40} {: <40} {: <40} {: <40} {: <40}'.format(
                                                                        f'Out(x_{atom}_{label}, "%11.5f")',
                                                                        f'Out(y_{atom}_{label}, "%11.5f")',
                                                                        f'Out(z_{atom}_{label}, "%11.5f")',
                                                                        f'Out(occ_{atom}_{label}, "%11.5f")',
                                                                        f'Out(beq_{atom}_{label}, "%11.5f")',
            )
        )
        
        fout.write('\n')

    fout.write('\t\tOut_String("\\n")\n')
    fout.write('#endif')
    fout.write('\n\n')




def read_cif(path):

    data = {'atoms': {}} # Initialise dictionary
    read = True # Initialise read toggle

    # Lists attributes to get out of the .CIF-file. This will correspond to what VESTA writes out, not necessarily what you will find in ICSD
    attrs = ['_chemical_name_common', '_cell', '_space_group_name_H-M_alt', '_space_group_IT_number']

    # Open file
    with open(path, 'r') as cif:
        line = cif.readline()

        # Read until encountering #End
        while read:
            
            # Break loop if #End is reached
            if not line or line.startswith('#End'):
                read = False
                break

            # Handle loops
            if line.lstrip().startswith('loop_'):
                loop = []
                line = cif.readline()

                # Only handle loops with attributes that starts with _atom. Other loops are incompatible with below code due to slight differences in formatting (e.g. lineshifts)
                while line.lstrip().startswith('_atom'):
                    loop.append(line) # add the attributes of the loop to a list (every keywords starting with _atom)
                    line = cif.readline() # Read next line

                
                # If there were any attributes that started with _atom - if this is empty is just means there were other attributes in this loop
                if loop:
                    # Read every line after the attribute listing has ended - until it encounters a new attribute, loop or #End tag
                    # FIXME WHat a horrible condition statement - need to fix this!
                    while line and not line.lstrip().startswith('_') and not line.lstrip().startswith('loop_') and not line.lstrip().startswith('#End') and not line=='\n':
                        # Initialise empty dictionary for a given atom if it has not already been created in another loop
                        if line.split()[0] not in data['atoms'].keys():
                            data["atoms"][line.split()[0]] = {}

                        # Add all the attribute / value pairs for the current loop
                        for i, attr in enumerate(loop):
                            data["atoms"][line.split()[0]][attr[:-1].lstrip()] = line.split()[i]
                        
                        # Read new line
                        line = cif.readline()
                
                # If loop list was empty, keep going to the next line
                else:
                    line = cif.readline()

            

            # Handle free-standing attributes
            elif any([line.lstrip().startswith(i) for i in attrs]):
            #line.startswith() or line.startswith('_symmetry_space_group_name') or line.startswith('_symmetry_Int_Tables'):
            #
                #
                attr, *value = line.split()

                value = ' '.join([str(i) for i in value])

                data[attr] = value
                line = cif.readline()

            else:
                line = cif.readline()   


    print(data.keys()) 
    return data
        



def update_topas_options(options):
    import nafuma.xrd as xrd
    basedir = os.path.dirname(xrd.__file__)

    with open(os.path.join(basedir, 'topas_default.json'), 'r') as default:
        topas_defaults = json.load(default)

    options['topas_options'] = aux.update_options(options=options['topas_options'], required_options=topas_defaults.keys(), default_options=topas_defaults)

    return options['topas_options']

def make_big_inp(data: dict, options={}):
    ''' Generates a big .INP-file with all filenames found in data["path"]. Uses a template .INP-file (which has to be generated manually from an initial refinement in TOPAS) and appends this to a large .INP-file
    while changing the filenames. '''


    # FIXME Strip headers from initial INP file before copying it. 

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

    with open(options['template'], 'r') as template:
        strip_headers(template)

    
    with open('tmp_template.inp', 'r') as template, open(options['output'], 'w', newline='\n') as output, open(runlist, 'w', newline='\n') as runlist:
    
        write_headers(output, options)

        template = template.read()
        
        for i, path in enumerate(data['path']):

            s = make_inp_entry(template=template, xdd=path, num=i, options=options)
            output.write(s)
            output.write('\n\n')
            
            runlist.write('#define \tUSE_'+f'{i}'.zfill(4) + '\n')

    os.remove('tmp_template.inp')


def strip_headers(fin):

    line = fin.readline()
    newlines = []

    while 'r_wp' not in line:
        if line[0] != '\'':
            line = fin.readline()            
        else:
            newlines.append(line)
            line = fin.readline()

    newlines.append(line)
    newlines = newlines + fin.readlines()

    with open('tmp_template.inp', 'w') as fout:
        for line in newlines:
            fout.write(line)

    







def write_headers(fout, options):
    # FIXME Could modify this to make sure that certain combinations of options is not written, such as both do_errors and bootstrap_errors. 

    headers = [
        "A_matrix_memory_allowed_in_Mbytes", 
        "approximate_A",
        "bootstrap_errors",
        "conserve_memory",
        "continue_after_convergence",
        "do_errors",
        "chi2_convergence_criteria", 
        "iters",
        "num_runs" ]


    for file in options['include']:
        fout.write(f'#include {file} \n')

    fout.write('\n')

    if options['save_results']:
        fout.write('#define output \n')
    
    fout.write('\n')

    for option, value in options['topas_options'].items():
        if value and option in headers:
            if isinstance(value, bool):
                fout.write(f'{option} \n')
            else:
                fout.write(f'{option} {value} \n')


def get_headers(inp, path):

    with open(inp, 'r') as inp:
        headers = ['index']

        line = inp.readline()

        while not path in line:
            line = inp.readline()

        # Jump down to lines
        line = inp.readline()
        line = inp.readline()

        while not '#endif' in line:
            if line.split():
                
                regx = r"\([\S]*"
                headers_line = re.findall(regx, line)

                for i, header in enumerate(headers_line):
                    header = header[1:-1]

                    if all(keyword in header for keyword in ['Get', '(', ')']):
                        header = header[4:-1]

                    headers_line[i] = header

                for header in headers_line:
                    if header != '"\\n"':
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

    # FIXME Since the big INP files now have the same filename for all iterations, we need to adjust the code to only get unique values from the get_paths function
    # FIXME get_headers() is also not working now. Needs to be adjusted to the new way of writing the Out-parameters    
    paths = get_paths(data['inp'])
    paths = aux.get_unique(paths)
    


    for path in paths:
        headers = get_headers(data['inp'], path)

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





def read_results(path):
    
    results = pd.read_csv(path, delim_whitespace=True, index_col=0)
    
    return results


