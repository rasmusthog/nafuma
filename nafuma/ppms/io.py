import pandas as pd
import numpy as np

import nafuma.auxillary as aux

def read_data(path, options={}):

    index = find_start(path)

    df = pd.read_csv(path, skiprows=index+1)
    mask = df.loc[df['Comment'].notna()]

    df = df[['Comment', 'Time Stamp (sec)', 'Temperature (K)', 'Magnetic Field (Oe)', 
    'DC Moment (emu)', 'DC Std. Err. (emu)', 'DC Quad. Moment (emu)', 
    'AC=1 DC=2 Locate=3', 'Max. Field (Oe)', 'Pressure (Torr)', 'Temp. Status (code)',
    ]]

    new_columns = ['Comment', 'Time', 'Temperature', 'Magnetic_Field', 
    'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 
    'Status', 'Max_Field', 'Pressure', 'Temperature_Status']

    df.columns = new_columns

    df[['Temperature', 'Magnetic_Field', 'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 'Max_Field', 'Pressure']] = df[['Temperature', 'Magnetic_Field', 'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 'Max_Field', 'Pressure']].astype(float)

    df = df.loc[df['DC_Std_Err'] < 0.001]

    if all([option in options.keys() for option in ['molar_mass', 'sample_mass']]):
        df = calculate_emu_per_mol_oe(df, options)
        df = calculate_bohr_magnetons(df, options)
        df = calculate_chi_inverse(df, options)


    dfs = []
    for i in range(1,len(mask.index)):
        dfs.append(df.iloc[mask.index[i-1]:mask.index[i]])

    return dfs


def read_hysteresis(path):

    index = find_start(path)

    df = pd.read_csv(path, skiprows=index+1)

    df = df[['Comment', 'Time Stamp (sec)', 'Temperature (K)', 'Magnetic Field (Oe)', 
    'DC Moment (emu)', 'DC Std. Err. (emu)', 'DC Quad. Moment (emu)', 
    'AC=1 DC=2 Locate=3', 'Max. Field (Oe)', 'Pressure (Torr)', 'Temp. Status (code)',
    ]]

    new_columns = ['Comment', 'Time', 'Temperature', 'Magnetic_Field', 
    'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 
    'Status', 'Max_Field', 'Pressure', 'Temperature_Status']

    df.columns = new_columns

    df[['Temperature', 'Magnetic_Field', 'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 'Max_Field', 'Pressure']] = df[['Temperature', 'Magnetic_Field', 'DC_Moment', 'DC_Std_Err', 'DC_Quad_Moment', 'Max_Field', 'Pressure']].astype(float)

    df = df.loc[df['DC_Std_Err'] < 0.001]

    return df


def find_start(path):

    with open(path, 'r') as f:

        i = 0
        line = f.readline()

        while '[Data]' not in line:
            line = f.readline()
            i += 1


            if i > 1000:
                break


    return i



def calculate_emu_per_mol_oe(df, options={}):

    m = options['sample_mass'] / 1000 # convert from mg to g
    n = m / options['molar_mass']


    df['DC_Moment_emu_per_mol'] = df['DC_Moment'] / n
    df['DC_Moment_emu_per_mol_oe'] = df['DC_Moment'] / (n * df['Magnetic_Field'])


    return df




def calculate_bohr_magnetons(df, options={}):


    default_options = {
        'units': 'cgs',
    }

    options = aux.update_options(options=options, default_options=default_options)

    if options['units'] == 'cgs':
        df['bohr_magnetons'] = df['DC_Moment_emu_per_mol'] * 1.07828E20 / 6.023E23    ## mu_B per emu divided by Avogadro's number

    return df
    


def calculate_chi_inverse(df, options={}):

    df['chi_inverse'] = 1/ df['DC_Moment_emu_per_mol']

    return df