import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

import warnings

import os

import linecache

import nafuma.auxillary as aux


def open_doscar(doscar='DOSCAR'):
    with open(doscar, 'r') as dos:
        lines = dos.readlines()
    
    return lines

def open_poscar(poscar='POSCAR'):
    with open(poscar, 'r') as pos:
        lines = pos.readlines()
        
    return lines

def open_outcar(outcar='OUTCAR'):
    with open(outcar, 'r') as out:
        lines = out.readlines()

    return lines

def get_number_of_atoms(doscar='DOSCAR'):
    lines = open_doscar(doscar)
    
    return int(lines[0].strip().split()[0])
    
def get_atoms(poscar):

    with open(poscar, 'r') as poscar:
        lines = poscar.readlines()

        atoms = lines[5].split()
        atom_num = lines[6].split()

        atom_num = [int(num) for num in atom_num]


    atoms_list = make_atoms_list(atoms, atom_num)

    atoms_info = {
        'specie': atoms, # list of specie in sequential order
        'number': atom_num, # list of number of per specie in sequential order
        'list': atoms_list # list of every individual atom in sequential order
    }

    return atoms_info
    
    
def get_fermi_level(doscar='DOSCAR', headerlines=6):
    lines = open_doscar(doscar)
    
    return float(lines[headerlines-1].strip().split()[3])


def get_valence_electron_count(outcar='OUTCAR', poscar='POSCAR'):
    lines = open_outcar(outcar)

    atoms, atoms_dict = get_atoms(poscar)
    n_atoms = len(atoms)

    for line in lines:
        line = line.strip()
        if line[0:4] == "ZVAL":
            valence_electrons = line.split()[-n_atoms:]
            break

    return valence_electrons

    


def read_coop(data={}, options={}):
    ''' Reads a COOPCAR.lobster file and prepares the DataFrame for further data handling.

    Input: 
    path: The path to the COOPCAR.lobster-file

    Output:
    coopcar: A DataFrame containing the COOP data
    interactions: A list of interactions'''


    required_options = ['collapse']

    default_options = {
        'collapse': False,
        'adjust': None,
    }


    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
     
    interactions = determine_interactions(data['coopcar'])

    coopcar = pd.read_csv(data['coopcar'], skiprows=3+len(interactions), header=None, delim_whitespace=True)

    spin_polarised = determine_spin_polarisation(coopcar, interactions)


    # Create list of column names 
    # If the run is spin polarised, add all up states first and then all down states 
    if spin_polarised:
        columns = ['Energy', 'avg_up', 'avg_int_up']
        for i in range(1, len(interactions)+1):
            columns += ['interaction{}_up'.format(i), 'interaction{}_int_up'.format(i)]
        columns += ['avg_down', 'avg_int_down']
        for i in range(1, len(interactions)+1):
            columns += ['interaction{}_down'.format(i), 'interaction{}_int_down'.format(i)]

    # Otherwise just 
    else:
        columns = ['Energy', 'avg', 'avg_int']
        for i in range(1, len(interactions)+1):
            columns += ['interaction_{}'.format(i), 'interaction_{}_int'.format(i)]




    coopcar.columns = columns

    if options['adjust']:
        coopcar['Energy'] = coopcar['Energy'] - options['adjust']


    if options['collapse']:
        
        columns_collapsed = ['Energy']

        for column in columns:
            if column.split('_')[0] not in columns_collapsed:
                columns_collapsed.append(''.join(column.split('_')[:-1]))
                columns_collapsed.append(''.join(column.split('_')[:-1])+'_int')


        for column in columns_collapsed:
            if column != 'Energy':
                coopcar[column] = coopcar[column+'_up'] + coopcar[column+'_down']

        columns.remove('Energy')
        coopcar.drop(columns, axis=1, inplace=True)
        coopcar.columns = columns_collapsed





    return coopcar, interactions


def read_cohp(path, flip_sign=False):
    ''' Reads a COHPCAR.lobster file and prepares the DataFrame for further data handling.

    Input: 
    path: The path to the COOPCAR.lobster-file
    flip_sign: Boolean value to determine whether all COHP-values should be flipped (that is, return -COHP)

    Output:
    coopcar: A DataFrame containing the COOP data
    interactions: A list of interactions'''

    interactions = determine_interactions(path)

    cohpcar = pd.read_csv(path, skiprows=3+len(interactions), header=None, delim_whitespace=True)

    spin_polarised = determine_spin_polarisation(cohpcar, interactions)


    # Create list of column names 
    # If the run is spin polarised, add all up states first and then all down states 
    if spin_polarised:
        columns = ['Energy', 'avg_up', 'avg_up_int']
        for i in range(1, len(interactions)+1):
            columns += ['interaction_{}_up'.format(i), 'interaction_{}_up_int'.format(i)]
        columns += ['avg_down', 'avg_down_int']
        for i in range(1, len(interactions)+1):
            columns += ['interaction_{}_down'.format(i), 'interaction_{}_up_int'.format(i)]

    # Otherwise just 
    else:
        columns = ['Energy', 'avg', 'avg_int']
        for i in range(1, len(interactions)+1):
            columns += ['interaction_{}'.format(i), 'interaction_{}_int'.format(i)]


    cohpcar.columns = columns

    if flip_sign:
        columns = columns[1:]

    for column in columns:
        cohpcar[column] = -cohpcar[column]

    return cohpcar, interactions


def determine_interactions(path):
    ''' Determines the number of interactions present in the COOPCAR.lobster file. 

    Input:
    path: The path to the COOPCAR.lobster-file

    Output:
    interactions: A list of strings with the interactions in the COOPCAR.lobster file'''


    with open(path, 'r') as coop:
        lines = coop.readlines()


    interactions = []
    for line in lines:
        if line[0:2] == 'No':
            interactions.append(line.split(':')[-1].split('(')[0])


    return interactions



def determine_spin_polarisation(coopcar, interactions):
    ''' Determines whether a COOPCAR.lobster file is spin polarised or not.

    Input:
    coopcar: A DataFrame containing the COOP data
    interactions: A list of interactions obtained from the determine_interactions function

    Output:
    spin_polarised: Boolean value to indicate whether or not the COOP data is spin polarised'''

    number_of_interactions = len(interactions)

    spin_polarised = True if coopcar.shape[1] == 4*number_of_interactions + 5 else False

    return spin_polarised




def read_dos(path, flip_down=True):

    dos_info = get_doscar_information(path)

    with open(path, 'r') as doscar:
        count = 0
        raw_data = []

        while count < dos_info["NEDOS"] + 6:
            if count >= 6:
                data_line = [float(x) for x in doscar.readline().split()]
                raw_data.append(data_line)

            else:
                doscar.readline()
            
            count += 1


    dos = pd.DataFrame(raw_data)

    if dos_info["spin_polarised"]:
        header_names = ['energy', 'total_up', 'total_down', 'total_integrated_up', 'total_integrated_down']
    else:
        header_names = ['energy', 'total', 'total_integrated']


    dos.columns = header_names

    if dos_info["spin_polarised"] and flip_down:
        dos["total_down"] = -dos["total_down"]
        dos["total_integrated_down"] = -dos["total_integrated_down"]

    return dos


def read_pdos(data: dict, options={}):
    ''' data-dictionary should be structured like this:
    data["path"] - dictionary with path to POSCAR/CONTCAR mapped to key 'poscar' and path to DOSCAR in key 'doscar'. '''

    default_options = {
        'sum_atoms': False,
        'sum_orbitals': False,
        'adjust': False,        # Manually adjust the energy scale if the automatic Fermi-level detection didn't work.
        'normalise': False,
    }

    options = aux.update_options(options=options, default_options=default_options)

    # Grab some metadata
    data['atoms'] = get_atoms(data['path']['poscar'])
    data['info'] = get_doscar_information(data['path']['doscar'])

    
    # Read the data
    with open(data['path']['doscar'], 'r') as f:
        doscar_raw = f.readlines()

    # Read basis functions - this will only yield a non-empty list if DOSCAR is generated by LOBSTER. It will only contain columns for the used basis set, so to be able to make proper column headers, this information is needed.
    data['basis_functions'] = get_basis_functions(doscar_raw)

    # If DOSCAR is from VASP and not LOBSTER, the basis functions will not be listed in DOSCAR (but all columns will be present)
    if not data['basis_functions']:
        data['basis_functions'] = [None for _ in data['atoms']['list']]

    
    data['pdos'] = []
    for i, (atom, basis_functions) in enumerate(zip(data['atoms']['list'], data['basis_functions'])):

        # Define line to start read-in: (1+i) to skip total DOS in the start, plus every previous PDOS, including the headerline. Then the initial 6 header lines in the start of the file.
        start = int( (1 + i) * ( data['info']['NEDOS'] + 1 ) + 6 )

        pdos_atom = []
        for j in range(0,data['info']['NEDOS']):
            pdos_atom.append(doscar_raw[start+j].split())

        # Convert the data to a DataFrame and convert to float
        pdos_atom = pd.DataFrame(pdos_atom)
        pdos_atom = pdos_atom.astype(float)


        # If DOSCAR was written by VASP, set the standard columns to the basis_functions variable
        # FIXME This does not allow for f-orbitals. Could do a check vs. shape of pdos_atom to determine
        if not basis_functions:
            basis_functions = ['s', 'p1', 'p2', 'p3', '2p1', '2p2', '2p3', 'd1', 'd2', 'd3', 'd4', 'd5']

    
        # Split the basis functions into up and down spin channels
        if data['info']['spin_polarised']:
            columns = ['Energy']
            for function in basis_functions:
                columns.append(function+'_up')
                columns.append(function+'_down')

        else:
            columns = ['Energy'] + basis_functions
        
        # Set the new column names
        pdos_atom.columns = columns

        if options['adjust']:
            pdos_atom['Energy'] -= options['adjust']

        if options['smooth']:
            pdos_atom = smoothing(pdos=pdos_atom, options=options)

        if options['normalise']:
            orbitals = [orbital for orbital in columns if 'Energy' not in orbital]
            for k, specie in enumerate(data['atoms']['specie']):
                if specie == atom:
                    index = k

            for orbital in orbitals:
                pdos_atom[orbital] = pdos_atom[orbital] / data['atoms']['number'][index]

        # Make total columns
        if data['info']['spin_polarised']:
            up_functions = [orbital for orbital in columns if '_up' in orbital]
            down_functions = [orbital for orbital in columns if '_down' in orbital]

            pdos_atom['total_up'] = pdos_atom[up_functions].sum(axis=1)
            pdos_atom['total_down'] = pdos_atom[down_functions].sum(axis=1)
            pdos_atom['total'] = pdos_atom[['total_up', 'total_down']].sum(axis=1)

            # Flip all the values for spin down
            down_functions.append('total_down')
            for orbital in down_functions:
                pdos_atom[orbital] = (-1)*pdos_atom[orbital]

        data['pdos'].append(pdos_atom)

    if options['sum_atoms']:
        data['pdos'] = sum_atoms(data, options)
    if options['sum_orbitals']:
        data['pdos'] = sum_orbitals(data, options)

    return data


def sum_atoms(data: dict, options={}):
    ''' Needs to have read data through read_pdos() '''

    default_options = {

    }

    options = aux.update_options(options=options, default_options=default_options)


    total_atoms = 0

    # Sort all the DataFrames into a dictionary
    sorted_dfs = {}
    for i, (specie, number) in enumerate(zip(data['atoms']['specie'], data['atoms']['number'])):
        
        sorted_dfs[specie] = []

        for j in range(number):
            index = j if i == 0 else j + total_atoms
            
            sorted_dfs[specie].append(data['pdos'][index])

        total_atoms += number

    # Sum all the DataFrames for each specie
    summed_dfs = {}
    for i, (specie, dfs) in enumerate(sorted_dfs.items()):
        for j, df in enumerate(dfs):
            if j > 0:
                dfs[0] += df

        dfs[0]['Energy'] = dfs[0]['Energy'] / data['atoms']['number'][i]
        summed_dfs[specie] = dfs[0]

    return summed_dfs


def sum_orbitals(data: dict, options={}):
    
    shells = [
        '1s', '2s', '3s', '4s', '5s', '6s',
        '2p', '3p', '4p', '5p', '6p',
        '3d', '4d', '5d', '6d',
        '4f', '5f', '6f'
    ]

    summed_dfs = []

    # If sum_atoms() is already called
    if isinstance(data['pdos'], dict):
        summed_dfs = {}

        for specie, specie_pdos in data['pdos'].items():
            orbitals = [orbital for orbital in specie_pdos.columns if 'Energy' not in orbital]
            new_pdos = pd.DataFrame(specie_pdos['Energy'])

            for shell in shells:
                if data['info']['spin_polarised']:
                    sub_columns_up = [orbital for orbital in orbitals if orbital.startswith(shell) and orbital.endswith('_up')]
                    sub_columns_down = [orbital for orbital in orbitals if orbital.startswith(shell) and orbital.endswith('_down')]
                    
                    new_pdos[shell+'_up'] = data['pdos'][specie][sub_columns_up].sum(axis=1)
                    new_pdos[shell+'_down'] = data['pdos'][specie][sub_columns_down].sum(axis=1)
                
                else:
                    sub_columns = [orbital for orbital in orbitals if orbital.startswith(shell)]
                    
                    new_pdos[shell] = data['pdos'][specie][sub_columns].sum(axis=1)
            new_pdos = new_pdos.loc[:, (new_pdos != 0).any(axis=0)]
            summed_dfs[specie] = new_pdos
                

    else:
        for pdos in data['pdos']:
            orbitals = [orbital for orbital in pdos.columns if 'Energy' not in orbital]
            new_pdos = pd.DataFrame(pdos['Energy'])
            for shell in shells:
                if data['info']['spin_polarised']:
                    sub_columns_up = [orbital for orbital in orbitals if orbital.startswith(shell) and orbital.endswith('_up')]
                    sub_columns_down = [orbital for orbital in orbitals if orbital.startswith(shell) and orbital.endswith('_down')]
                    new_pdos[shell+'_up'] = pdos[sub_columns_up].sum(axis=1)
                    new_pdos[shell+'_down'] = pdos[sub_columns_down].sum(axis=1)
                else:
                    sub_columns = [orbital for orbital in orbitals if orbital.startswith(shell)]
                    new_pdos[shell] = pdos[sub_columns].sum(axis=1)
            new_pdos = new_pdos.loc[:, (new_pdos != 0).any(axis=0)]
            summed_dfs.append(new_pdos)
    
    
    return summed_dfs


def smoothing(pdos: pd.DataFrame, options={}):
    ' Smoothes the data using the Savitzky-Golay filter. This is the only algorithm at this moment.  '

    default_options = {
        'smooth_window_length': 3, # Determines the window length of smoothing that the savgol-filter uses for smoothing
        'smooth_polyorder': 2, # Determines the order of the polynomial used in the smoothing algorithm
        'smooth_algorithm': 'savgol', # At the present, only Savitzky-Golay filter is implemented. Add Gaussian and Boxcar later.
    }

    options = aux.update_options(options=options, default_options=default_options)

    # Initialise new DataFrame with correct x-values
    pdos_smooth = pd.DataFrame(pdos['Energy'])

    orbitals = [orbital for orbital in pdos.columns if 'Energy' not in orbital]

    if options['smooth_algorithm'] == 'savgol':
        
        for orbital in orbitals:
            pdos_smooth.insert(1, orbital, savgol_filter(pdos[orbital], options['smooth_window_length'], options['smooth_polyorder']))

    return pdos_smooth



def get_basis_functions(doscar):

    basis_functions = []
    for line in doscar:
        if 'Z=' in line:
            basis_functions.append(line.split(';')[-1].split())


    return basis_functions



def get_doscar_information(path):
    ''' Reads information from the DOSCAR'''

    kind = 'LOBSTER' if 'LOBSTER' in linecache.getline(path, 5) else 'VASP'
    dos_info_raw = linecache.getline(path, 6).split()
    spin_polarised = True if len(linecache.getline(path, 7).split()) == 5 else False


    dos_info = {
        'ENMIN': float(dos_info_raw[0]), 
        'ENMAX': float(dos_info_raw[1]), 
        'NEDOS': int(dos_info_raw[2]), 
        'EFERMI': float(dos_info_raw[3]), 
        'spin_polarised': spin_polarised, 
        'kind': kind
    }

    return dos_info


def make_atoms_list(atoms, number_of_atoms):

    atoms_list = []
    for atom, num in zip(atoms, number_of_atoms):
        atoms_list = atoms_list + [atom for _ in range(num)]

    return atoms_list
        


def read_pdos_legacy(data: dict, options={}):

    required_options = ['flip_down', 'sum_atoms', 'sum_orbitals', 'collapse_spin', 'adjust', 'manual_adjust', 'normalise', 'normalise_unit_cell', 'normalisation_factor']

    default_options = {
        'flip_down': True,
        'sum_atoms': False,
        'sum_orbitals': False,
        'collapse_spin': False,
        'adjust': False,
        'manual_adjust': None,
        'normalise': False,
        'normalise_unit_cell': False,
        'normalisation_factor': None
    }
    
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    # Get information about the DOSCAR-file (ENMIN, ENMAX, NEDOS, EFERMI and spin_polarised)
    dos_info = get_doscar_information(data['doscar'])

    # Get information from the POSCAR-file (or CONTCAR) (which species, number of atoms per specie and a dictionary matching the two)
    species, atom_num, atoms_dict = get_atoms(data['poscar'])


    # Open DOSCAR and read all lines
    with open(data['doscar'], 'r') as file:
        doscar_data = file.readlines()

    # Make list of all individual atoms
    atoms = []
    for specie in species:
        for i in range(0,atoms_dict[specie]):
            atoms.append(specie)


    # Loop through all the atoms and make a DataFrame for each atom.

    pdos_full = []
    for ind, atom in enumerate(atoms):

        # Define line to start reading
        start = int((1 + ind) * dos_info["NEDOS"] + 6 + (ind * 1))
        
        add_d_orbitals = False
        add_p_orbitals = False


        
        
        # Check if d-orbitals are included (if DOSCAR comes from LOBSTER), otherwise they will have to be added to make the DataFrames the same size    
        if dos_info['kind'] == 'LOBSTER':
            # Extract basis sets from each atom
            basis_sets = doscar_data[start].split(';')[-1]

            if not 'd' in basis_sets:
                add_d_orbitals = True

            if basis_sets.count('p') != 6:
                add_p_orbitals = True
                

        
        # Start reading lines into a list and convert to DataFrame and cast all entries as float
        pdos_atom = []
        for i in range(1,dos_info["NEDOS"]+1):
            pdos_atom.append(doscar_data[start+i].split())


        pdos_atom = pd.DataFrame(pdos_atom)
        pdos_atom = pdos_atom.astype(float)




        # Give the columns names and add second set of p-orbitals, d-orbitals, or both if they don't exist in the data
        if add_p_orbitals:
            if dos_info['spin_polarised']:
                pdos_atom.columns = ['Energy', 's_u', 's_d', 'p1_u', 'p1_d', 'p2_u', 'p2_d', 'p3_u', 'p3_d']
                pdos_atom[['2p1_u', '2p1_d', '2p2_u', '2p2_d', '2p3_u', '2p3_d']] = 0

                if add_d_orbitals:
                    pdos_atom[['d1_u', 'd1_d', 'd2_u', 'd2_d', 'd3_u', 'd3_d', 'd4_u', 'd4_d', 'd5_u', 'd5_d']] = 0
                    pdos_atom['tot_u'] = pdos_atom[['s_u', 'p1_u', 'p2_u', 'p3_u', '2p1_u', '2p2_u', '2p3_u', 'd1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                    pdos_atom['tot_d'] = pdos_atom[['s_d', 'p1_d', 'p2_d', 'p3_d', '2p1_d', '2p2_d', '2p3_d', 'd1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)

                else:
                    pdos_atom['tot_u'] = pdos_atom[['s_u', 'p1_u', 'p2_u', 'p3_u', '2p1_u', '2p2_u', '2p3_u', 'd1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                    pdos_atom['tot_d'] = pdos_atom[['s_d', 'p1_d', 'p2_d', 'p3_d', '2p1_d', '2p2_d', '2p3_d', 'd1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)



        elif add_d_orbitals:
            if dos_info["spin_polarised"]:
                pdos_atom.columns = ['Energy', 's_u', 's_d', 'p1_u', 'p1_d', 'p2_u', 'p2_d', 'p3_u', 'p3_d']
                pdos_atom[['d1_u', 'd1_d', 'd2_u', 'd2_d', 'd3_u', 'd3_d', 'd4_u', 'd4_d', 'd5_u', 'd5_d']] = 0
                pdos_atom['tot_u'] = pdos_atom[['s_u', 'p1_u', 'p2_u', 'p3_u', '2p1_u', '2p2_u', '2p3_u', 'd1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                pdos_atom['tot_d'] = pdos_atom[['s_d', 'p1_d', 'p2_d', 'p3_d', '2p1_d', '2p2_d', '2p3_d', 'd1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)


            else:
                pdos_atom.columns = ['Energy', 's', 'p1', 'p2', 'p3']
                pdos_atom[['d1', 'd2', 'd3', 'd4', 'd5']] = 0
                pdos_atom['tot'] = pdos_atom[['s', 'p1', 'p2', 'p3', '2p1', '2p2', '2p3', 'd1', 'd2', 'd3', 'd4', 'd5']].sum(axis=1)

        else:
            if dos_info["spin_polarised"]:
                if dos_info['kind'] == 'LOBSTER':
                    pdos_atom.columns = ['Energy', 's_u', 's_d', 'p1_u', 'p1_d', 'p2_u', 'p2_d', 'p3_u', 'p3_d', '2p1_u', '2p1_d', '2p2_u', '2p2_d', '2p3_u', '2p3_d', 'd1_u', 'd1_d', 'd2_u', 'd2_d', 'd3_u', 'd3_d', 'd4_u', 'd4_d', 'd5_u', 'd5_d']
                    pdos_atom['tot_u'] = pdos_atom[['s_u', 'p1_u', 'p2_u', 'p3_u', '2p1_u', '2p2_u', '2p3_u','d1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                    pdos_atom['tot_d'] = pdos_atom[['s_d', 'p1_d', 'p2_d', 'p3_d', '2p1_d', '2p2_d', '2p3_d', 'd1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)

                elif dos_info['kind'] == 'VASP':
                    pdos_atom.columns = ['Energy', 's_u', 's_d', 'p1_u', 'p1_d', 'p2_u', 'p2_d', 'p3_u', 'p3_d', 'd1_u', 'd1_d', 'd2_u', 'd2_d', 'd3_u', 'd3_d', 'd4_u', 'd4_d', 'd5_u', 'd5_d']
                    pdos_atom['tot_u'] = pdos_atom[['s_u', 'p1_u', 'p2_u', 'p3_u', 'd1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                    pdos_atom['tot_d'] = pdos_atom[['s_d', 'p1_d', 'p2_d', 'p3_d', 'd1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)

            else:
                pdos_atom.columns = ['Energy', 's', 'p1', 'p2', 'p3', '2p1', '2p2', '2p3', 'd1', 'd2', 'd3', 'd4', 'd5']
                pdos_atom['tot'] = pdos_atom[['s', 'p1', 'p2', 'p3', '2p1', '2p2', '2p3', 'd1', 'd2', 'd3', 'd4', 'd5']].sum(axis=1)



        # Change the sign of the spin down columns if flip_down is True
        if options['flip_down'] and dos_info['spin_polarised'] and not options['collapse_spin']:
            down = [orbital for orbital in pdos_atom.columns if '_d' in orbital]

            for orbital in down:
                pdos_atom[orbital] = -pdos_atom[orbital]


        if options['manual_adjust']:
            pdos_atom["Energy"] = pdos_atom["Energy"] - options['manual_adjust']



        pdos_full.append(pdos_atom)


    # If sum_atoms is True, all DataFrames of a given specie will be added together
    if options['sum_atoms']:

        # Initalise a new emtpy list that will become our pdos_full
        pdos_full_sum_atoms = []


        start = 0
        # Loop through each specie so that each specie will get exactly one DataFrame
        for specie in species:
            # Initialise with first DataFrame of the specie
            atom_pdos_summed = pdos_full[start]

            # Loop over DOSes and add if there's a match for specie
            for ind, pdos in enumerate(pdos_full):
                if atoms[ind] == specie and ind != start:
                    atom_pdos_summed = atom_pdos_summed + pdos


            # Divide the Energy by the number of DataFrames added to get back to the original value of the Energy
            atom_pdos_summed["Energy"] = atom_pdos_summed["Energy"] / atoms_dict[specie]


            if options['normalise']:

                if dos_info["spin_polarised"]:
                    columns = atom_pdos_summed.columns
                    #columns = ['s_u', 's_d', 'p1_u', 'p1_d', 'p2_u', 'p2_d', 'p3_u', 'p3_d', '2p1_u', '2p1_d', '2p2_u', '2p2_d', '2p3_u', '2p3_d', 'd1_u', 'd1_d', 'd2_u', 'd2_d', 'd3_u', 'd3_d', 'd4_u', 'd4_d', 'd5_u', 'd5_d', 'tot_u', 'tot_d']
                    for column in columns:
                        atom_pdos_summed[column] = atom_pdos_summed[column] / atoms_dict[specie]


            # Append the new DataFrame for a given specie to the list
            pdos_full_sum_atoms.append(atom_pdos_summed)

            start += atoms_dict[specie]


        # Rename the list
        pdos_full = pdos_full_sum_atoms


    # If collapse_spin is True for a spin polarised DOSCAR, the up and down channels of each orbital will be added together. 
    if options['collapse_spin'] and dos_info['spin_polarised']:
        
        pdos_full_spin_collapsed = []
        for pdos in pdos_full:
            temp_pdos = pd.DataFrame()

            temp_pdos["Energy"] = pdos["Energy"]

            temp_pdos['s'] =  pdos[['s_u', 's_d']].sum(axis=1)

            temp_pdos['p1'] = pdos[['p1_u', 'p1_d']].sum(axis=1)
            temp_pdos['p2'] = pdos[['p2_u', 'p2_d']].sum(axis=1)
            temp_pdos['p3'] = pdos[['p3_u', 'p3_d']].sum(axis=1)

            temp_pdos['2p1'] = pdos[['2p1_u', '2p1_d']].sum(axis=1)
            temp_pdos['2p2'] = pdos[['2p2_u', '2p2_d']].sum(axis=1)
            temp_pdos['2p3'] = pdos[['2p3_u', '2p3_d']].sum(axis=1)

            temp_pdos['d1'] = pdos[['d1_u', 'd1_d']].sum(axis=1)
            temp_pdos['d2'] = pdos[['d2_u', 'd2_d']].sum(axis=1)
            temp_pdos['d3'] = pdos[['d3_u', 'd3_d']].sum(axis=1)
            temp_pdos['d4'] = pdos[['d4_u', 'd4_d']].sum(axis=1)
            temp_pdos['d5'] = pdos[['d5_u', 'd5_d']].sum(axis=1)

            temp_pdos['tot'] = pdos[['tot_u', 'tot_d']].sum(axis=1)

            pdos_full_spin_collapsed.append(temp_pdos)


        pdos_full = pdos_full_spin_collapsed
        dos_info['spin_polarised'] = False



    # If sum_orbitals is True, all columns belonging to a particular set of orbitals will be added together.
    if options['sum_orbitals']:
        pdos_full_sum_orbitals = []
        
        for pdos in pdos_full:
            temp_pdos = pd.DataFrame()

            temp_pdos["Energy"] = pdos["Energy"]

            if dos_info['spin_polarised']:
                temp_pdos['s_u'] = pdos['s_u']
                temp_pdos['s_d'] = pdos['s_d']
                
                temp_pdos['p_u'] = pdos[['p1_u', 'p2_u', 'p3_u']].sum(axis=1)
                temp_pdos['p_d'] = pdos[['p1_d', 'p2_d', 'p3_d']].sum(axis=1)

                if len(pdos_full[0].columns) == 25:
                    temp_pdos['2p_u'] = pdos[['2p1_u', '2p2_u', '2p3_u']].sum(axis=1)
                    temp_pdos['2p_d'] = pdos[['2p1_d', '2p2_d', '2p3_d']].sum(axis=1)

                temp_pdos['d_u'] = pdos[['d1_u', 'd2_u', 'd3_u', 'd4_u', 'd5_u']].sum(axis=1)
                temp_pdos['d_d'] = pdos[['d1_d', 'd2_d', 'd3_d', 'd4_d', 'd5_d']].sum(axis=1)

                temp_pdos['tot_u'] = pdos['tot_u']
                temp_pdos['tot_d'] = pdos['tot_d']

            else:
                temp_pdos['s'] = pdos['s']
                temp_pdos['p'] = pdos[['p1', 'p2', 'p3']].sum(axis=1)
                temp_pdos['2p'] = pdos[['2p1', '2p2', '2p3']].sum(axis=1)
                temp_pdos['d'] = pdos[['d1', 'd2', 'd3', 'd4', 'd5']].sum(axis=1)
                temp_pdos['tot'] = pdos['tot']

            pdos_full_sum_orbitals.append(temp_pdos)

        pdos_full = pdos_full_sum_orbitals




    return pdos_full, dos_info





#def get_bader_charges(poscar='POSCAR', acf='ACF.dat'):



def write_poscar(data: dict, options={}):

    required_options = ['path', 'overwrite', 'scale']

    default_options = {
        'path': './POSCAR.vasp',
        'overwrite': False,
        'scale': 1.0
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if os.path.isfile(options['path']) and not options['overwrite']:
        warnings.warn(f"File {options['path']} already exists, and overwrite disabled. Quitting.")
        return None


    atom_nums = count_atoms(data)

    with open(options['path'], 'w') as fout:
        
        # Write first line
        for atom in data['atoms']:
            if atom_nums[atom] > 1:
                fout.write(f'{atom}{atom_nums[atom]}')
            else:
                fout.write(f'{atom}')

        fout.write('\t - Automatically generated by the NAFUMA Python package.\n')

        # Write second line
        fout.write(str(options['scale'])+'\n')

        # Write lattice parameters
        # FIXME Now only writes cells without any angles
        fout.write('\t{:<09} \t {:<09} \t {:<09}\n'.format(
                                                        str(data['lattice_parameters'][0]),
                                                        str(0.0),
                                                        str(0.0),
                                                    )
                                                )

        fout.write('\t{:<09} \t {:<09} \t {:<09}\n'.format(
                                                        str(0.0),
                                                        str(data['lattice_parameters'][1]),
                                                        str(0.0),
                                                    )
                                                )

        fout.write('\t{:<09} \t {:<09} \t {:<09}\n'.format(
                                                        str(0.0),
                                                        str(0.0),
                                                        str(data['lattice_parameters'][2]),
                                                    )
                                                )


        # Write atoms
        for atom in data['atoms']:
            fout.write(f'\t{atom}')
        fout.write('\n')

        for atom in data['atoms']:
            fout.write(f'\t{atom_nums[atom]}')
        fout.write('\n')

        fout.write('Direct\n')



    
        for atom in data['atoms']:
            for label, coords in data['coordinates'].items():
                if atom in label:
                    fout.write('\t{:<09} \t {:<09} \t {:<09}\n'.format(
                                                                    coords[0],
                                                                    coords[1],
                                                                    coords[2]
                                                                )
                                                            )


    

def apply_transformation(data, rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], translation=[0,0,0]):

    if not isinstance(rotation, np.ndarray):
        rotation = np.array(rotation)

    if not rotation.shape == (3,3):
        print("NOOOO!!!!")

    for label in data['coordinates'].keys():
        data['coordinates'][label] = rotation.dot(data['coordinates'][label])
        data['coordinates'][label] = translate_back(data['coordinates'][label])
        data['coordinates'][label] = data['coordinates'][label] + translation

    return data


def translate_back(coords):

    for i, coord in enumerate(coords):
        if coord < 0:
            while coords[i] < 0:
                coords[i] = coords[i]+1

        elif coord >= 1:
            while coords[i] >= 1:
                coords[i] = coords[i]-1

    return coords


def count_atoms(data):
    atom_nums = {}
    for atom in data['atoms']:
        atom_nums[atom] = 0

    for label in data['coordinates'].keys():
        for atom in data['atoms']:
            if atom in label:
                atom_nums[atom] += 1


    return atom_nums

def append_data(data, new_data):

    if not new_data:
        return data

    if not data:
        data = {
            'atoms': new_data['atoms'],
            'coordinates': {},
            'lattice_parameters': new_data['lattice_parameters']
        }

    atom_num = count_atoms(data)

    new_coords = {}

    for label, coords in data['coordinates'].items():
        new_coords[label] = coords


    extend_unit_cell = [0,0,0]
    for label, coords in new_data['coordinates'].items():
        atom = ''.join(filter(str.isalpha, label))
        number = int(''.join(filter(str.isnumeric, label)))
        new_number = number + atom_num[atom]
        new_label =  atom+str(new_number)

        new_coords[new_label] = coords

        for i, coord in enumerate(coords):
            if coord > 1 and np.floor(coord) >= extend_unit_cell[i]:
                extend_unit_cell[i] = np.floor(coord)

    data['coordinates'] = new_coords

    return data


def make_supercell(data, supercell):

    for i, param in enumerate(data['lattice_parameters']):
        data['lattice_parameters'][i] = supercell[i] * param

        if supercell[i] > 0:

            for label, coords in data['coordinates'].items():
                data['coordinates'][label][i] = (1/supercell[i])*data['coordinates'][label][i]



    
def copy_data(data):

    new_data = {}

    new_data['atoms'] = list(data['atoms'])
    new_data['coordinates'] = data['coordinates'].copy()
    new_data['lattice_parameters'] = list(data['lattice_parameters'])

    return new_data


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))