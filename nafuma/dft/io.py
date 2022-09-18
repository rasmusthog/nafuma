import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    atoms_dict = {}

    for ind, atom in enumerate(atoms):
        atoms_dict[atom] = atom_num[ind]

    return atoms, atom_num, atoms_dict
    
    
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
        'collapse': False
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


def get_doscar_information(path):
    ''' Reads information from the DOSCAR'''

    kind = 'LOBSTER' if 'LOBSTER' in linecache.getline(path, 5) else 'VASP'
    dos_info_raw = linecache.getline(path, 6).split()
    spin_polarised = True if len(linecache.getline(path, 7).split()) == 5 else False


    dos_info = {'ENMIN': float(dos_info_raw[0]), 'ENMAX': float(dos_info_raw[1]), 'NEDOS': int(dos_info_raw[2]), 'EFERMI': float(dos_info_raw[3]), 'spin_polarised': spin_polarised, 'kind': kind}

    return dos_info


#def get_bader_charges(poscar='POSCAR', acf='ACF.dat'):

