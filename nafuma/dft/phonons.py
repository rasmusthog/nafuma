import re
from this import d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import subprocess
import os
import shutil

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import importlib
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from cycler import cycler
import itertools


def get_atoms(path='.'):

    poscar = os.path.join(path, 'POSCAR')

    with open(poscar, 'r') as poscar:
        lines = poscar.readlines()

        atoms = lines[5].split()
        atom_num = lines[6].split()


        atom_num = [int(num) for num in atom_num]

    return atoms, atom_num

def get_dimensions(path='.'):

    poscar = os.path.join(path, 'POSCAR')
    sposcar = os.path.join(path, 'SPOSCAR')


    with open(poscar, 'r') as poscar:

        lines_pos = poscar.readlines()


    with open(sposcar, 'r') as sposcar:

        lines_spos = sposcar.readlines()



    a_p, b_p, c_p = lines_pos[2].split(), lines_pos[3].split(), lines_pos[4].split()
    a_s, b_s, c_s = lines_spos[2].split(), lines_spos[3].split(), lines_spos[4].split()

    lattice_params_poscar, lattice_params_sposcar = [a_p, b_p, c_p], [a_s, b_s, c_s]



    poscar_new = []
    sposcar_new = []

    for lp_p, lp_s in zip(lattice_params_poscar, lattice_params_sposcar):
        lp_p = np.sqrt(float(lp_p[0])**2 + float(lp_p[1])**2 + float(lp_p[2])**2)
        lp_s = np.sqrt(float(lp_s[0])**2 + float(lp_s[1])**2 + float(lp_s[2])**2)

        poscar_new.append(lp_p)
        sposcar_new.append(lp_s)


    dim = [int(lp_s/lp_p) for lp_s, lp_p in zip(sposcar_new, poscar_new)]

    return dim


def read_band(band_dir):
    ''' Reads a band file as written by the function write_phonon_bands() into a pandas DataFrame and returns this. Contains two columns: k-points (the "distance" output in the band.yaml-file by phonopy) and frequencies.
    
    Input:
    band_dir: the path to the band-file.
    
    Output:
    band: pandas DataFrame containing frequencies of the band along the k-point path specified in the phonopy calculation'''
    
    
    # Read the band into a pandas DataFrame
    band = pd.read_csv(band_dir, delim_whitespace=True, header=None, names=['kpt', 'frequency'])
    
    return band


def read_kpoints(kpoints_dir):
    ''' Reads a VASP KPOINTS-file in line mode. Returns two lists: special_points_coords, containing the coordinates of the special points in k-space and special_points_labels, the names of these special points.
    Requires a KPOINTS-file that is in line mode with special points indicated with a "!". 
    
    Input:
    kpoints_dir: the path to the KPOINTS-file
    
    Output:
    special_points_coords: List of 3D coordinates of the k-space special points
    special_points_labels: List of names of the k-space special points'''
    
    
    # Open the KPOINTS-file and read it line by line, appending each line with a "!" to the special_points list.
    special_points = []
    
    with open(kpoints_dir) as kpoints:
        lines = kpoints.readlines()
        
        for line in lines:
            if '!' in line:
                special_points.append(line)
    
    
    # Go through the special points to separate them into the coordinate and the label for each special point into special_points_coords and special_points_labels respectively.
    special_points_coords = []
    special_points_labels = []
    
    for special_point in special_points:
        if len(special_point.split()) == 5:
            special_points_coords.append(special_point.split()[0:3])
            special_points_labels.append(special_point.split()[-1])
                
    
    return special_points_coords, special_points_labels





def get_kpoints_ticks(band):
    ''' Finds the coordinates for the special points in the 1D-projection given by phonopy (the parameter 'distance' in the band.yaml file). This is to determine the placement of labels and vertical lines in the bandstructure plot.
    
    Input:
    band: the path to a band_XX.dat file. Should not matter which one is passed here.
    
    Output:
    kpts_ticks: A list of coordinates corresponding to the special points.'''
    
    band = np.genfromtxt(band)
    
    kpts_ticks = []
    
    # Append the first point
    kpts_ticks.append(0.)
    
    # Go through all data points - where the x-value repeats, a k-point tick is appended to the list
    for j in np.arange(np.shape(band)[0]-1):
        if band[j,0]==band[j+1,0]:
            kpts_ticks.append(band[j,0])
        
    # Append the last point
    kpts_ticks.append(max(band[:,0]))
    
    
    return kpts_ticks
                
def get_kpoints_labels(special_points_labels):
    ''' Takes the raw special point labels from read_kpoints() and writes them in a way to be used in the bandstructure plots. 
    Where there is a discontinuity in the path, the label is separated with a |. 
    
    Input:
    special_points_labels: A list of special points as directly read from the KPOINTS-file by read_kpoints()
    
    Ouput:
    labels: A list of labels suitable to pass as x-ticks during plotting of the bandstructure plots.'''
    
    
    # Loop through the raw special points labels list following a set of rules, to extract the labels suitable for plotting
    labels = []
    
    for ind, label in enumerate(special_points_labels):
        
        # Add the first label as this will be a separate special point 
        if ind == 0:
            label = '${}$'.format(label) if (label[0] == '\\') else label
            labels.append(label)
            
        # Add the last label, as this will also be a separate special point (or will it? Must change this if that is not always the case)
        elif ind == len(special_points_labels)-1:
            label = '${}$'.format(label) if (label[0] == '\\') else label
            labels.append(label)
            
        # Skip every second entry, as they will repeat due to the way the KPOINTS-file is constructed
        elif ind%2 != 0:
            continue
        
        # Add label if it's continuous (i.e. if the current and previous points are the same), add "previous|current" if discontinuous (i.e. if they are not the same)
        else:
            if label == special_points_labels[ind-1]:
                label = '${}$'.format(label) if (label[0] == '\\') else label # If the special point has a greek letter, such as the gamma point, makes sure that the label is enclosed in $ to be rendered correctly.
                labels.append(label)
                
            else:
                label = '${}$'.format(label) if (label[0] == '\\') else label # If the special point has a greek letter, such as the gamma point, makes sure that the label is enclosed in $ to be rendered correctly.
                previous_label = special_points_labels[ind-1]
                previous_label = '${}$'.format(previous_label) if previous_label[0] == '\\' else previous_label # If the special point has a greek letter, such as the gamma point, makes sure that the label is enclosed in $ to be rendered correctly.
                
                labels.append("{}|{}".format(previous_label, label))
                       
    
    return labels
    
    
def read_phonon_dos(dos_path):
    ''' Reads the phonon density of states from a total_dos.dat file as written by phonopy. This file will be generated by the function calculate_phonon_dos() as well as this calls phonopy to calculate the density of states.
    
    Input:
    dos_path: the path to the total_dos.dat file. Must include the filename
    
    Output:
    df: pandas DataFrame containing the contents of the total_dos.dat file. Two columns, "Frequency" and "DOS". '''
    
    df = pd.read_csv(dos_path, header=None, skiprows=1, delim_whitespace=True)
    df.columns = ['Frequency', 'DOS']
    
    return df


def read_phonon_pdos(path, normalise=False, poscar=None):
    ''' Reads the phonon density of states from a total_dos.dat file as written by phonopy. This file will be generated by the function calculate_phonon_dos() as well as this calls phonopy to calculate the density of states.
    
    Input:
    dos_path: the path to the total_dos.dat file. Must include the filename
    
    Output:
    df: pandas DataFrame containing the contents of the total_dos.dat file. Two columns, "Frequency" and "DOS". '''
    
    df = pd.read_csv(path, index_col=0)

    if normalise and poscar:
        atoms, atom_num = get_atoms(poscar)


        for atom, num in zip(atoms, atom_num):
            df[atom] = df[atom] / num

    
    return df
    
    
def write_phonopy_band_path(special_points_coords):
    ''' Writes the band path used by phonopy to calculate the bandstructure from the raw information as extracted by read_kpoints(). 
    
    Input:
    special_points_coords: list of coordinates for the special points as read by read_kpoints(), that reads a VASP KPOINTS.bands file.
    
    Output:
    phonopy_band_path: '''
    
    coords = []
    
    for ind, coord in enumerate(special_points_coords):
        
        # Add the first label
        if ind == 0:
            coord = "{} {} {} ".format(coord[0], coord[1], coord[2])
            coords.append(coord)
            
        # Add the last label
        elif ind == len(special_points_coords)-1:
            coord = "{} {} {}".format(coord[0], coord[1], coord[2])
            coords.append(coord)
            
        # Skip every second entry
        elif ind%2 != 0:
            continue
        
        # Add label if it's continuous, add "previous|current" if discontinuous
        else:
            if coord == special_points_coords[ind-1]:
                coord = "{} {} {} ".format(coord[0], coord[1], coord[2])
                coords.append(coord)
                
            else:
                first_coord = "{} {} {}".format(special_points_coords[ind-1][0], special_points_coords[ind-1][1], special_points_coords[ind-1][2])
                second_coord = "{} {} {} ".format(coord[0], coord[1], coord[2])
                coords.append("{}, {} ".format(first_coord, second_coord))


    phonopy_band_path = ''
    
    for coord in coords:
        phonopy_band_path = phonopy_band_path + coord
                
    return phonopy_band_path


def write_mesh_conf(atoms, dim, mesh, dos_range=None, pdos=False, atom_num=None, tmax=None):
    
    atom_str = 'ATOM_NAME = '
    for atom in atoms:
        atom_str += atom + " "
        
    dim_str = 'DIM = '
    for d in dim:
        dim_str += str(d) + " "
        
    mesh_str = 'MP = '
    for m in mesh:
        mesh_str += str(m) + " "
        
    dos_str = 'DOS_RANGE = '
    for d in dos_range:
        dos_str += str(d) + " "

    if tmax:
        tmax_str = f'TMAX = {tmax}'

    if pdos:
        pdos_str = 'PDOS ='
    
        atoms_sum = 0
        for ind, atom in enumerate(atoms):
            for i in range(1,atom_num[ind]+1):
                pdos_str += " {}".format(i+atoms_sum)
                
            # Add comma after numbers unless it's the last entry
            if ind != len(atom_num)-1:
                pdos_str += ','
                
            atoms_sum = atoms_sum + atom_num[ind]
        
        
    with open('mesh.conf', 'w') as conf:
        
        conf.write(atom_str + '\n' + dim_str + '\n' + mesh_str)
            
        if tmax:
            conf.write('\n' + tmax_str)
        
        if dos_range:
            conf.write('\n' + dos_str)
            
        if pdos:
            conf.write('\n' + pdos_str)
            
        conf.write('\n' + "WRITE_MESH = .FALSE.")
        

def write_band_conf(atoms, dim, mesh, band, band_points=None):
    
    atom_str = 'ATOM_NAME = '
    for atom in atoms:
        atom_str += atom + " "
        
    dim_str = 'DIM = '
    for d in dim:
        dim_str += str(d) + " "
        
    mesh_str = 'MP = '
    for m in mesh:
        mesh_str += str(m) + " "
        
    band_str = 'BAND = ' + band
    
    if band_points:
        band_points_str = "BAND_POINTS = " + band_points
    
    
    with open('band.conf', 'w') as conf:
        
        if not band_points:
            conf.write(atom_str + '\n' + dim_str + '\n' + mesh_str + '\n' + band_str)
            
        else:
            conf.write(atom_str + '\n' + dim_str + '\n' + mesh_str + '\n' + band_str + '\n' + band_points_str)
            
            
            
def calculate_phonon_dos(path, atoms, dim, mesh, dos_range=None):
    
    cwd = os.getcwd()
    os.chdir(path)

    write_mesh_conf(atoms, dim, mesh, dos_range=dos_range)
    
    subprocess.call('phonopy -ps mesh.conf >> phonopy_output.dat', shell=True)


    # Make folder and move output in there
    os.mkdir('total_dos')
    shutil.move('total_dos.pdf', 'total_dos/total_dos.pdf')
    shutil.move('total_dos.dat', 'total_dos/total_dos.dat')
    shutil.move('mesh.conf', 'total_dos/mesh.conf')
    shutil.move('phonopy.yaml', 'total_dos/phonopy.yaml')
    shutil.move('phonopy_output.dat', 'total_dos/phonopy_output.dat')

    os.chdir(cwd)





def calculate_phonon_pdos(path, atoms, dim, mesh, dos_range=None, atom_num=None, order=None):
    ''' Calculate the projected phonon DOS. Calls function to write mesh.conf file, and then cleans up the output by summing all the individual contributions per atom to the same species.'''
    
    cwd = os.getcwd()
    os.chdir(path)

    write_mesh_conf(atoms, dim, mesh, dos_range=dos_range, pdos=True, atom_num=atom_num)
    
    subprocess.call('phonopy -ps mesh.conf >> phonopy_output.dat', shell=True)


    df = pd.read_csv('projected_dos.dat', delim_whitespace=True, skiprows=1, header=None, dtype=float)


    # Loop over the columns and add according to "atoms" and "atom_num" lists

    atoms_sum = 0
    for atom, num in zip(atoms, atom_num):
        df[atom] = df[1+atoms_sum]
    
        for i in range(2+atoms_sum, atoms_sum+num+1):
            df[atom] = df[atom] + df[i]
        
        
        
        atoms_sum += num
        
       

    # Remove all other columns, and rename the first column to "Frequency" 
    df.drop(df.iloc[:, 1:atoms_sum+1], inplace = True, axis = 1)
    df.rename(columns = {0: "Frequency"}, inplace=True)


    # If a list is passed to order, this will change the order of the atoms:

    if order:
        df_temp = pd.DataFrame()
        df_temp["Frequency"] = df["Frequency"]

        for atom in order:
            df_temp[atom] = df[atom]

        df = df_temp


    # Save the cleaned up DataFrame to file.
    df.to_csv('projected_dos_clean.dat')
    


    # Make folder and move output in there
    os.mkdir('projected_dos')
    shutil.move('partial_dos.pdf', 'projected_dos/partial_dos.pdf')
    shutil.move('projected_dos.dat', 'projected_dos/projected_dos.dat')
    shutil.move('projected_dos_clean.dat', 'projected_dos/projected_dos_clean.dat')
    shutil.move('mesh.conf', 'projected_dos/mesh.conf')
    shutil.move('phonopy.yaml', 'projected_dos/phonopy.yaml')
    shutil.move('phonopy_output.dat', 'projected_dos/phonopy_output.dat')



    os.chdir(cwd)



def calculate_thermal_properties(path, atoms, dim, mesh, dos_range=None, tmax=None):
    
    cwd = os.getcwd()
    os.chdir(path)

    write_mesh_conf(atoms, dim, mesh, dos_range=dos_range, tmax=tmax)
    
    subprocess.call('phonopy -t mesh.conf >> phonopy_output.dat', shell=True)

    with open('phonopy_output.dat', 'r') as f:
        lines = f.readlines()


    data = []
    for ind, line in enumerate(lines):

        if line.split():
            if "#" in line.split()[0]:
                j = 1
                while lines[ind+j].split():
                    data.append(lines[ind+j].split())
                    j += 1



    df = pd.DataFrame(data)
    df.columns = ['T', 'F', 'S', 'C_v', 'E']

    df.to_csv('thermal_properties.dat')




    #Make folder and move output in there
    os.mkdir('thermal_properties')
    shutil.move('thermal_properties.yaml', 'thermal_properties/thermal_properties.yaml')
    shutil.move('thermal_properties.dat', 'thermal_properties/thermal_properties.dat')
    shutil.move('mesh.conf', 'thermal_properties/mesh.conf')
    shutil.move('phonopy.yaml', 'thermal_properties/phonopy.yaml')
    shutil.move('phonopy_output.dat', 'thermal_properties/phonopy_output.dat')
    
    
    os.chdir(cwd)
            
def calculate_phonon_bandstructure(path, atoms, dim, mesh, kpoints='KPOINTS.bands', band_points=None):
    
    cwd = os.getcwd()
    os.chdir(path)


    kpoints_coords, kpoints_labels = read_kpoints(kpoints)
    
    band = write_phonopy_band_path(kpoints_coords)
    
    write_band_conf(atoms, dim, mesh, band, band_points=band_points)
    
    subprocess.call('phonopy band.conf >> phonopy_output.dat', shell=True)
    
    write_phonon_bands()


    os.mkdir('dispersion_relation')

    shutil.move('band.conf', 'dispersion_relation/band.conf')
    shutil.move('band.yaml', 'dispersion_relation/band.yaml')
    shutil.move('bands', 'dispersion_relation/bands')
    shutil.move('mesh.yaml', 'dispersion_relation/mesh.yaml')
    shutil.move('phonopy.yaml', 'dispersion_relation/phonopy.yaml')
    shutil.move('phonopy_output.dat', 'dispersion_relation/phonopy_output.dat')
    
    os.chdir(cwd)

def write_phonon_bands(band='band.yaml'):
    
    with open(band, 'r') as f:
        lines = f.readlines()
    
    
    kpoints = []
    frequencies = []
    
    for line in lines:
        if 'distance' in line:
            kpoints.append(line.split()[-1])
            
        if 'frequency' in line:
            frequencies.append(line.split()[-1])

            
    number_of_kpoints = len(kpoints)
    number_of_bands = len(frequencies) / number_of_kpoints
    
    if not os.path.isdir('bands'):
        os.mkdir('bands')
        
    os.chdir('bands')
    
    for i in range(int(number_of_bands)):
        
        with open('band_{}.dat'.format(i+1), 'w') as b:
            for ind, kpoint in enumerate(kpoints):
                if ind == len(kpoints)-1:
                    b.write("{} {}".format(kpoint, frequencies[ind*int(number_of_bands)+i]))
                else:
                    b.write("{} {}\n".format(kpoint, frequencies[ind*int(number_of_bands)+i]))
                
                
    os.chdir('../')
    

def plot_phonon_dos(dos_path='total_dos.dat', options={}):
    
    
    required_options = ['xlim', 'ylim', 'flip_xy', 'colours', 'palettes', 'rc_params', 'format_params']


    default_options = {
        'xlim': None, # x-limits
        'ylim': None, # y-limits
        'flip_xy': False,  # Whether to flip what is plotted on the x- and y-axes respectively. Default is False and plots frequency along x-axis and density of states along y-axis.
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'format_params': {},
        'rc_params': {}       
    }


    options = update_options(options=options, required_options=required_options, default_options=default_options)
    
    fig, ax = prepare_plot(options=options)

    dos = read_phonon_dos(dos_path=dos_path)
    
    if not options['xlim']:
        options['xlim'] = [dos["Frequency"].min(), dos["Frequency"].max()]
    
    if not options['ylim']:
        options['ylim'] = [dos["DOS"].min(), dos["DOS"].max()*1.1]
    
    
    if not options['colours']:
        colours = generate_colours(palette=options['palette'])
    else:
        colours = itertools.cycle(options['colours'])
    
    if options['flip_xy']:
        dos.plot(x='DOS', y='Frequency', ax=ax, color=colours[0])
        
    else:
        dos.plot(x='Frequency', y='DOS', ax=ax, color=colours[0])
    
    
    options['plot_kind'] = 'DOS'
    fig, ax = prettify_dos_plot(fig=fig, ax=ax, options=options)
    
    ax.get_legend().remove()


    return fig, ax




def prepare_plot(options={}):

    rc_params = options['rc_params']
    format_params = options['format_params']
    
    required_options = ['single_column_width', 'double_column_width', 'column_type', 'width_ratio', 'aspect_ratio', 'compress_width', 'compress_height', 'upscaling_factor', 'dpi']

    default_options = {
    'single_column_width': 8.3,
    'double_column_width': 17.1,
    'column_type': 'single',
    'width_ratio': '1:1',
    'aspect_ratio': '1:1',
    'compress_width': 1,
    'compress_height': 1,
    'upscaling_factor': 1.0,
    'dpi': 600, 
    }
    
    options = update_options(format_params, required_options, default_options)


    # Reset run commands
    plt.rcdefaults()
    
    # Update run commands if any is passed (will pass an empty dictionary if not passed)
    update_rc_params(rc_params)
    
    width = determine_width(options)
    height = determine_height(options, width)
    width, height = scale_figure(options=options, width=width, height=height)
    
    fig, ax = plt.subplots(figsize=(width, height), dpi=options['dpi'])
    
    return fig, ax


def plot_phonon_pdos(path='projected_dos_clean.dat', options={}):
    

    required_options = ['xlim', 'ylim', 'flip_xy', 'colours', 'palettes', 'normalise', 'poscar', 'atoms', 'rc_params', 'format_params']


    default_options = {
        'xlim': None, # x-limits
        'ylim': None, # y-limits
        'flip_xy': False,  # Whether to flip what is plotted on the x- and y-axes respectively. Default is False and plots frequency along x-axis and density of states along y-axis.
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'normalise': False,
        'poscar': None,
        'atoms': [],
        'format_params': {},
        'rc_params': {}
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    dos = read_phonon_pdos(path=path, normalise=options['normalise'], poscar=options['poscar'])
    
    fig, ax = prepare_plot(options=options)
    
    if not options['xlim']:
        options['xlim'] = [dos["Frequency"].min(), dos["Frequency"].max()]
    
    if not options['ylim'] and options['atoms']:
        ymin = 0
        ymax = 0


        for atom in options['atoms']:
            if dos[atom].min() < ymin:
                ymin = dos[atom].min()

            if dos[atom].max() > ymax:
                ymax = dos[atom].max()

        options['ylim'] = [ymin, ymax*1.1]
    
    
    if not options['colours']:
        colours = generate_colours(palette=options['palette'])
    else:
        colours = itertools.cycle(options['colours'])

    for ind, atom in enumerate(options['atoms']):
    
        if options['flip_xy']:
            dos.plot(x=atom, y='Frequency', ax=ax, color=next(colours))
        
        else:
            dos.plot(x='Frequency', y=atom, ax=ax, color=next(colours))
    
    
    options['plot_kind'] = 'PDOS'
    prettify_dos_plot(fig=fig, ax=ax, options=options)
    
    return fig, ax

    

    
def plot_phonon_bandstructure(band_folder='bands', kpoints='KPOINTS.bands', options={}, title=None, xlim=None, ylim=None, pad_bottom=None, scale=1, square=True, width=None, height=None, dpi=None, rotation=None, xpad=None, ypad=None):
    
    
    # Get the special points labels
    kpoint_coords, kpoint_labels = read_kpoints(kpoints)
    kpoint_labels = get_kpoints_labels(kpoint_labels)
    
    
    
    # Get current folder and change into the folder containing bands
    cwd = os.getcwd()
    os.chdir(band_folder)
    
    band_paths = [band for band in os.listdir() if os.path.isfile(band) and band[0:4] == 'band']
    
    # Get the location of the special points along the x-axis
    kpoint_ticks = get_kpoints_ticks(band_paths[0])
    
    bands = []
    for band_path in band_paths:
        bands.append(read_band(band_path))
    


    fig, ax = prepare_plot(options=options)
    
    mod = importlib.import_module("palettable.colorbrewer.%s" % 'qualitative')
    colour = getattr(mod, 'Dark2_3').mpl_colors[0]
    
    kpt_min = None
    kpt_max = None
    freq_min = None
    freq_max = None
    
    for band in bands:
        if kpt_min == None or band["kpt"].min() < kpt_min:
            kpt_min = band["kpt"].min()
        if kpt_max == None or band["kpt"].max() > kpt_max:
            kpt_max = band["kpt"].max()

        if freq_min == None or band["frequency"].min() < freq_min:
            freq_min = band["frequency"].min()
        if freq_max == None or band["frequency"].max() > freq_max:
            freq_max = band["frequency"].max()
            
        band.plot('kpt', 'frequency', ax=ax, color=colour)


    if not xlim:
        xlim = [kpt_min, kpt_max]
    if not ylim:
        ylim = [freq_min-freq_max*0.1, freq_max+freq_max*0.1]
    
    ax.get_legend().remove()
    
    prettify_plot(fig=fig, ax=ax, special_points_labels=kpoint_labels, special_points_coords=kpoint_ticks, xlim=xlim, ylim=ylim, title=title, pad_bottom=pad_bottom, scale=scale, rotation=rotation, xpad=xpad, ypad=ypad)
    
    os.chdir(cwd)


def prepare_plot_old(width=None, height=None, square=True, dpi=None, colour_cycle=('qualitative', 'Dark2_8'), temperatureunit='K', energyunit='eV f.u.$^{-1}$', scale=1):
    
    linewidth = 3*scale
    axeswidth = 3*scale
    
    plt.rc('lines', linewidth=linewidth)
    plt.rc('axes', linewidth=axeswidth)
    
    if square:
        if not width:
            width = 20
    
        if not height:
            height = width
        
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width, height), facecolor='w', dpi=dpi)
    
    return fig, ax


def prettify_plot(fig, ax, frequencyunit='THz', special_points_coords=None, special_points_labels=None, xlim=None, ylim=None, title=None, pad_bottom=None, scale=1, rotation=None, xpad=None, ypad=None):
    

    # Set sizes of ticks, labes etc.
    ticksize = 30*scale
    labelsize = 30*scale
    legendsize = 30*scale
    titlesize = 30*scale
    
    linewidth = 3*scale
    axeswidth = 3*scale
    majorticklength = 20*scale
    minorticklength = 10*scale
    

    # Set labels on x- and y-axes
    if ypad:
        ax.set_ylabel('Frequency [{}]'.format(frequencyunit), size=labelsize, labelpad=ypad)
    else:
        ax.set_ylabel('Frequency [{}]'.format(frequencyunit), size=labelsize)        


    ax.set_xlabel('')
    

    
    
    ax.tick_params(axis='y', direction='in', which='major', right=True, length=10, width=0.5)
    ax.tick_params(axis='y', direction='in', which='minor', right=True, length=5, width=0.5)
    
    ax.tick_params(axis='x', direction='in', which='major', bottom=False)
    
   
    
        
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(2.5))
                                      
    
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    
    
    
    # Set tick parameters
    if special_points_coords:
        for coord in special_points_coords:
            plt.axvline(coord, color='black', linestyle='--', linewidth=0.5)
    
    
    
    plt.xticks(ticks=special_points_coords, labels=special_points_labels, rotation=rotation)

    
    
    if xlim:
        plt.xlim(xlim)
    
    if ylim:
        plt.ylim(ylim)
        
        
    if title:
        ax.set_title(title, size=40)
        
    if pad_bottom is not None:
            bigax = fig.add_subplot(111) 
            bigax.set_facecolor([1,1,1,0])
            bigax.spines['top'].set_visible(False)
            bigax.spines['bottom'].set_visible(True)
            bigax.spines['left'].set_visible(False)
            bigax.spines['right'].set_visible(False)
            bigax.tick_params(labelcolor='w', color='w', direction='in', top=False, bottom=True, left=False, right=False, labelleft=False, pad=pad_bottom)

    if xpad:
        ax.tick_params(axis='x', pad=xpad)
    if ypad:
        ax.tick_params(axis='y', pad=ypad)
    
    return fig, ax


def prettify_dos_plot(fig, ax, options, frequencyunit='THz', dosunit='a.u.', xlim=None, ylim=None, title=None, hide_ylabels=False, flip_xy=False, pad_bottom=None, scale=1, pdos=False, colours=None, atoms=None, xpad=None, ypad=None):
    

    required_options = ['plot_kind', 'flip_xy', 'hide_x_labels', 'hide_y_labels',  'xlabel', 'ylabel', 'xunit', 'yunit', 'xlim', 'ylim', 'x_tick_locators', 'y_tick_locators', 'hide_x_ticks', 'hide_y_ticks', 'hide_x_ticklabels', 'hide_y_ticklabels',
                        'colours', 'palettes',  'title', 'legend', 'legend_position', 'subplots_adjust', 'text']

    default_options = {
        'plot_kind': 'DOS', # DOS or PDOS
        'flip_xy': False,
        'hide_x_labels': False, # Whether x labels should be hidden
        'hide_x_ticklabels': False,
        'hide_x_ticks': False,
        'hide_y_labels': False, # whether y labels should be hidden
        'hide_y_ticklabels': False,
        'hide_y_ticks': False,
        'xlabel': 'Frequency',
        'ylabel': 'DOS',
        'xunit': r'THz', # The unit of the x-values in the curve plot
        'yunit': r'a.u.', # The unit of the y-values in the curve and bar plots
        'xlim': None,
        'ylim': None,
        'x_tick_locators': [5, 2.5], # Major and minor tick locators
        'y_tick_locators': [10, 5],
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'title': None,
        'legend': True,
        'legend_position': ['upper center', (0.20, 0.90)], # the position of the legend passed as arguments to loc and bbox_to_anchor respectively 
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9],
        'text': None
    }


    if 'plot_kind' in options.keys():
        if 'ylabel' not in options.keys():
            if options['plot_kind'] == 'DOS':
                options['ylabel'] = 'Density of states'
            elif options['plot_kind'] == 'PDOS':
                options['ylabel'] = 'PDOS'


    options = update_options(options=options, required_options=required_options, default_options=default_options)

    
    if options['flip_xy']:

        # Switch all the x- and y-specific values
        options = swap_values(dict=options, key1='xlim', key2='ylim')
        options = swap_values(dict=options, key1='xunit', key2='yunit')
        options = swap_values(dict=options, key1='xlabel', key2='ylabel')
        options = swap_values(dict=options, key1='x_tick_locators', key2='y_tick_locators')
        options = swap_values(dict=options, key1='hide_x_labels', key2='hide_y_labels')

    # Set labels on x- and y-axes
    if not options['hide_y_labels']:
        ax.set_ylabel(f'{options["ylabel"]} [{options["yunit"]}]')
    else:
        ax.set_ylabel('')
        


    if not options['hide_x_labels']:
        ax.set_xlabel(f'{options["xlabel"]} [{options["xunit"]}]')
    else:
        ax.set_xlabel('')
        

    # Hide x- and y- ticklabels
    if options['hide_y_ticklabels']:
        ax.tick_params(axis='y', direction='in', which='both', labelleft=False, labelright=False)
    if options['hide_x_ticklabels']:
        ax.tick_params(axis='x', direction='in', which='both', labelbottom=False, labeltop=False)


    # Hide x- and y-ticks:
    if options['hide_y_ticks']:
        ax.tick_params(axis='y', direction='in', which='both', left=False, right=False)
    if options['hide_x_ticks']:
        ax.tick_params(axis='x', direction='in', which='both', bottom=False, top=False)



    # Set multiple locators
    ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
    ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))

    ax.xaxis.set_major_locator(MultipleLocator(options['x_tick_locators'][0]))
    ax.xaxis.set_minor_locator(MultipleLocator(options['x_tick_locators'][1]))

           
    # Set title
    if options['title']:
        ax.set_title(options['title'])


    # Generate colours
    if not options['colours']:
        colours = generate_colours(palette=options['palette'])
    else:
        colours = itertools.cycle(options['colours'])
        


    # Create legend

    if ax.get_legend():
                ax.get_legend().remove()

            
    if options['legend']:
        if options['plot_kind'] == 'PDOS' and options['atoms']:
        
            # Create legend
            patches = []
            for atom in options['atoms']:
                patches.append(mpatches.Patch(color=next(colours), label=atom))
        
            fig.legend(handles=patches, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1], frameon=False)

        

    # Adjust where the axes start within the figure. Default value is 10% in from the left and bottom edges. Used to make room for the plot within the figure size (to avoid using bbox_inches='tight' in the savefig-command, as this screws with plot dimensions)
    plt.subplots_adjust(left=options['subplots_adjust'][0], bottom=options['subplots_adjust'][1], right=options['subplots_adjust'][2], top=options['subplots_adjust'][3])


    # If limits for x- and y-axes is passed, sets these.
    if options['xlim'] is not None:
        ax.set_xlim(options['xlim'])

    if options['ylim'] is not None:
        ax.set_ylim(options['ylim'])


    # Add custom text
    if options['text']:
        plt.text(x=options['text'][1][0], y=options['text'][1][1], s=options['text'][0])
    
    return fig, ax





def read_thermal_properties(path, number_of_formula_units=None, convert=True):
    
    kJ = 6.2415064799632E+21
    Na = 6.0221415E+23
    
    thermal_properties = pd.read_csv(path, skiprows=1, index_col=0)
    thermal_properties.columns = ['T', 'F', 'S', 'Cv', 'E']

    
    if convert:
        thermal_properties.F = thermal_properties.F / Na * kJ
        thermal_properties.S = thermal_properties.S / Na * kJ
        thermal_properties.Cv = thermal_properties.Cv / Na * kJ
        thermal_properties.E = thermal_properties.E / Na * kJ
    
    if number_of_formula_units:
        thermal_properties.F = thermal_properties.F / number_of_formula_units
        thermal_properties.S = thermal_properties.S / number_of_formula_units
        thermal_properties.Cv = thermal_properties.Cv / number_of_formula_units
        thermal_properties.E = thermal_properties.E / number_of_formula_units
    
    
    
    return thermal_properties



def plot_thermal_properties(path, number_of_formula_units=None, convert=True):
    
    thermal_properties = read_thermal_properties(path=path, number_of_formula_units=number_of_formula_units, convert=convert)
    
    thermal_properties.plot(x='T', y=['F', 'S', 'Cv', 'E'])





def get_adjusted_energies(paths, equilibrium_energies, options={}):

    
    required_options = ['plot_kind', 'reference', 'number_of_formula_units', 'xlim', 'ylim', 'flip_xy', 'colours', 'palettes', 'normalise', 'poscar', 'atoms', 'rc_params', 'format_params']


    default_options = {
        'plot_kind': 'absolute',
        'reference': 0,
        'number_of_formula_units': None,
        'xlim': None, # x-limits
        'ylim': None, # y-limits
        'flip_xy': False,  # Whether to flip what is plotted on the x- and y-axes respectively. Default is False and plots frequency along x-axis and density of states along y-axis.
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'normalise': False,
        'poscar': None,
        'atoms': [],
        'format_params': {},
        'rc_params': {}
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    dfs = []
    
    if not options['number_of_formula_units']:
        options['number_of_formula_units'] = [None for i in range(len(paths))]
    
    for ind, path in enumerate(paths):
        df = read_thermal_properties(path, options['number_of_formula_units'][ind])
        dfs.append(df)


    for ind, df in enumerate(dfs):
        df["adjusted_energy"] = equilibrium_energies[ind] + df["F"]
        
    
    if options['plot_kind'] == 'difference':
        for ind, df in enumerate(dfs):
            df["reference_energy"] = dfs[options['reference']]["adjusted_energy"]
            df["difference_energy"] = df["adjusted_energy"] - df["reference_energy"]


    if options['plot_kind'] == 'relative':
        for ind, df in enumerate(dfs):
            df["reference_energy"] = dfs[options['reference']]["adjusted_energy"].iloc[0]
            df["relative_energy"] = df["adjusted_energy"] - df["reference_energy"]


    return dfs



def find_low_energy_structures_at_extremas(dfs):
    
    energy_low_T = -1
    low_T_ind = -1
    energy_high_T = -1
    high_T_ind = -1
    
    for ind, df in enumerate(dfs):
        if low_T_ind == -1:
            low_T_ind = ind
            energy_low_T = df['adjusted_energy'].loc[df['T'] == df['T'].min()].values[0]
        
        elif df['adjusted_energy'].loc[df['T'] == df['T'].min()].values[0] < energy_low_T:
            low_T_ind = ind
            energy_low_T = df['adjusted_energy'].loc[df['T'] == df['T'].min()].values[0]

        if high_T_ind == -1:
            high_T_ind = ind
            energy_high_T = df['adjusted_energy'].loc[df['T'] == df['T'].max()].values[0]
        
        elif df['adjusted_energy'].loc[df['T'] == df['T'].max()].values[0] < energy_high_T:
            high_T_ind = ind
            energy_high_T = df['adjusted_energy'].loc[df['T'] == df['T'].max()].values[0]



    return [low_T_ind, high_T_ind]

def find_intersection(dfs, ind1, ind2):

    intersection = -1

    for T in dfs[0]['T']:

        if dfs[ind2]['adjusted_energy'].loc[dfs[ind2]['T'] == T].values[0] < dfs[ind1]['adjusted_energy'].loc[dfs[ind1]['T'] == T].values[0]:
            intersection = T
            break

    
    return intersection





def plot_adjusted_energies(paths, equilibrium_energies, options={}):
    
    
    ''' This function plots the adjusted total energies of a set of structures given a set of thermal properties calculated using phonopy. 
    
    paths: List of paths (strings) to the .csv-files with thermal properties.
    equilibrium_energies: List of equilibrium energies (floats) of pristine calculations
    labels: List of labels (strings) to be shown in the plot
    mode: Whether to plot as a difference plot ("difference_plot") or absolute units ("absolute"). Defaults to absolute
    difference_reference: Index of which structure should serve as the reference. Defaults to 0. 
    number_of_formula_units: List of number of formula units per unit cell (int, float) to scale the data properly. Defaults to None, meaning to scaling.
    width: Width of the plot. Defaults to None, meaning standard width is used.
    width: Height of the plot. Defaults to None, meaning standard height is used.
    dpi: Dots per inch. Defaults to None, meaning standard dpi is used.
    colour_cycle: Tuple with type of colour scheme from the colorbrewer: http://jiffyclub.github.io/palettable/colorbrewer/
    temperatureunit: The unit to plot the temperature in. Only K implemented so far.
    energyunit: The unit to plot the energy in. Only eV per f.u. impleneted so far.
    inset: Whether or not there should be an inset. This is not very well implemented, and may cause issues. Defaults to False.
    inset_lims: The x-limits of the inset. Defaults to None, meaning it will just try to figure it out itself.
    '''

    required_options = ['plot_kind', 'reference', 'number_of_formula_units', 'labels', 'xlim', 'ylim', 'colours', 'palettes', 'linestyles', 'rc_params', 'format_params', 'inset_xlim', 'inset_ylim', 'draw_intersection_main', 'draw_intersection_inset', 'intersection_indices', 'intersection_lw']


    default_options = {
        'plot_kind': 'absolute',
        'reference': 0,
        'number_of_formula_units': None,
        'labels': None,
        'xlim': None, # x-limits
        'ylim': None, # y-limits
        'inset_xlim': None,
        'inset_ylim': None,
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'linestyles': ['solid', 'dotted', 'dashed'],
        'format_params': {},
        'rc_params': {},
        'draw_intersection_main': False,
        'draw_intersection_inset': False,
        'intersection_indices': None,
        'intersection_lw': None,
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    energy_dfs = get_adjusted_energies(paths=paths, equilibrium_energies=equilibrium_energies, options=options)

    fig, ax = prepare_plot(options=options)
    
    
    if not options['labels']:
        options['labels'] = ['_' for i in range(len(paths))]


    if not options['colours']:
        colours = generate_colours(palettes=options['palettes'])
    else:
        colours = itertools.cycle(options['colours'])
        

    linestyles = itertools.cycle(options['linestyles'])
    
    for df in energy_dfs:
        if options['plot_kind'] == 'difference':
            df.plot(x='T', y='difference_energy', ax=ax, ls=next(linestyles), c=next(colours))
        elif options['plot_kind'] == 'relative':
            df.plot(x='T', y='relative_energy', ax=ax, ls=next(linestyles), c=next(colours))
        elif options['plot_kind'] == 'absolute':
            df.plot(x='T', y='adjusted_energy', ax=ax, ls=next(linestyles), c=next(colours))
        
        ax.set_xlim([int(df["T"].min()), int(df["T"].max())])

    


    
    fig, ax = prettify_thermal_plot(fig=fig, ax=ax, options=options)
    
       

    if options['inset_xlim']:
        inset_ax = prepare_inset_axes(ax, options)

        if not options['colours']:
            colours = generate_colours(palettes=options['palettes'])
        else:
            colours = itertools.cycle(options['colours'])
        

        linestyles = itertools.cycle(options['linestyles'])
      
   
        for df in energy_dfs:
            if options['plot_kind'] =='absolute':
                y = 'adjusted_energy'
            elif options['plot_kind'] == 'relative':
                y = 'relative_energy'
            elif options['plot_kind'] == 'difference':
                y = 'difference_energy'

            df.loc[(df["T"] >= options['inset_xlim'][0]) & (df["T"] <= options['inset_xlim'][1])].plot(x='T', y=y, ax=inset_ax, ls=next(linestyles), c=next(colours))
            inset_ax.set_xlim([options['inset_xlim'][0], options['inset_xlim'][1]])
            
            if options['inset_ylim']:
                inset_ax.set_ylim([options['inset_ylim'][0], options['inset_ylim'][1]])
        
        inset_ax.get_legend().remove()
        inset_ax.set_xlabel('')


    if options['draw_intersection_main'] or options['draw_intersection_inset']:

        if not options['intersection_indices']:
            options['intersection_indices'] = find_low_energy_structures_at_extremas(energy_dfs)

        if not options['intersection_lw']:
            options['intersection_lw'] = plt.rcParams['lines.linewidth']

        intersection = find_intersection(energy_dfs, options['intersection_indices'][0], options['intersection_indices'][1])

        if options['draw_intersection_main']:
            ax.axvline(x=intersection, ls='dashed', c='black', lw=options['intersection_lw'])
        if options['draw_intersection_inset']:
            inset_ax.axvline(x=intersection, ls='dashed', c='black', lw=options['intersection_lw'])


    return fig, ax
        
    



def prepare_thermal_plot(width=None, height=None, dpi=None, colour_cycle=('qualitative', 'Dark2_8'), temperatureunit='K', energyunit='eV f.u.$^{-1}$', scale=1):
    
    linewidth = 3*scale
    axeswidth = 3*scale
    
    plt.rc('lines', linewidth=linewidth)
    plt.rc('axes', linewidth=axeswidth)
    
    if not width:
        width = 20
    
    if not height:
        height = width
        
    
    fig = plt.figure(figsize=(width, height), facecolor='w', dpi=dpi)
    ax = plt.gca()
    
    # Set colour cycle
    mod = importlib.import_module("palettable.colorbrewer.%s" % colour_cycle[0])
    colors = getattr(mod, colour_cycle[1]).mpl_colors
    ax.set_prop_cycle(cycler('color', colors))
    
    return fig, ax



def prettify_thermal_plot(fig, ax, options):
    
    required_options = ['plot_kind', 'hide_x_labels', 'hide_y_labels',  'rotation_x_ticks', 'rotation_y_ticks', 'xlabel', 'ylabel', 'xunit', 'yunit', 'xlim', 'ylim', 'x_tick_locators', 'y_tick_locators', 'hide_x_ticks', 'hide_y_ticks', 'hide_x_ticklabels', 'hide_y_ticklabels',
                        'colours', 'palettes',  'title', 'legend', 'legend_position', 'subplots_adjust', 'text']

    default_options = {
        'plot_kind': 'absolute', # absolute, relative, difference
        'hide_x_labels': False, # Whether x labels should be hidden
        'hide_x_ticklabels': False,
        'hide_x_ticks': False,
        'rotation_x_ticks': 0,
        'hide_y_labels': False, # whether y labels should be hidden
        'hide_y_ticklabels': False,
        'hide_y_ticks': False,
        'rotation_y_ticks': 0,
        'xlabel': 'Temperature',
        'ylabel': 'Energy',
        'xunit': r'K', # The unit of the x-values in the curve plot
        'yunit': r'eV f.u.$^{-1}$', # The unit of the y-values in the curve and bar plots
        'xlim': None,
        'ylim': None,
        'x_tick_locators': [100, 50], # Major and minor tick locators
        'y_tick_locators': [10, 5],
        'labels': None,
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'title': None,
        'legend': True,
        'legend_position': ['upper center', (0.20, 0.90)], # the position of the legend passed as arguments to loc and bbox_to_anchor respectively 
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9],
        'text': None
    }


    if 'plot_kind' in options.keys():
        if 'ylabel' not in options.keys():
            if options['plot_kind'] == 'absolute':
                options['ylabel'] = 'Energy'
            elif options['plot_kind'] == 'relative':
                options['ylabel'] = 'Relative energy'
            elif options['plot_kind'] == 'difference':
                options['ylabel'] = 'Energy difference'

        if 'y_tick_locators' not in options.keys():
            if options['plot_kind'] == 'absolute' or options['plot_kind'] == 'relative':
                options['y_tick_locators'] = [1, 0.5]
            elif options['plot_kind'] == 'difference':
                options['y_tick_locators'] = [0.1, 0.05]


    options = update_options(options=options, required_options=required_options, default_options=default_options)

    # Set labels on x- and y-axes
    if not options['hide_y_labels']:
        ax.set_ylabel(f'{options["ylabel"]} [{options["yunit"]}]')
    else:
        ax.set_ylabel('')
        
    if not options['hide_x_labels']:
        ax.set_xlabel(f'{options["xlabel"]} [{options["xunit"]}]')
    else:
        ax.set_xlabel('')


    # Set multiple locators
    ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
    ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))

    ax.xaxis.set_major_locator(MultipleLocator(options['x_tick_locators'][0]))
    ax.xaxis.set_minor_locator(MultipleLocator(options['x_tick_locators'][1]))
        
    # Hide x- and y- ticklabels
    if options['hide_y_ticklabels']:
        ax.tick_params(axis='y', direction='in', which='both', labelleft=False, labelright=False)
    else:
        plt.xticks(rotation=options['rotation_x_ticks'])
        #ax.set_xticklabels(ax.get_xticks(), rotation = options['rotation_x_ticks'])

    if options['hide_x_ticklabels']:
        ax.tick_params(axis='x', direction='in', which='both', labelbottom=False, labeltop=False)
    else:
        pass
        #ax.set_yticklabels(ax.get_yticks(), rotation = options['rotation_y_ticks'])


    # Hide x- and y-ticks:
    if options['hide_y_ticks']:
        ax.tick_params(axis='y', direction='in', which='both', left=False, right=False)
    if options['hide_x_ticks']:
        ax.tick_params(axis='x', direction='in', which='both', bottom=False, top=False)





           
    # Set title
    if options['title']:
        ax.set_title(options['title'])



        

    # Create legend

    if ax.get_legend():
        ax.get_legend().remove()

    
    if options['legend']:
        

        # Make palette and linestyles from original parameters
        if not options['colours']:
            colours = generate_colours(palettes=options['palettes'])
        else:
            colours = itertools.cycle(options['colours'])
        

        linestyles = itertools.cycle(options['linestyles'])
        
        # Create legend
        custom_lines = []
        active_labels = []

        for label in options['labels']:


            # Discard next linestyle and colour if label is _
            if label == '_':
                _ = next(colours)
                _ = next(linestyles)

            else:
                custom_lines.append(Line2D([0], [0], color=next(colours), ls=next(linestyles)))
                active_labels.append(label)

    

        ax.legend(custom_lines, active_labels, frameon=False, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1])
        #fig.legend(handles=patches, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1], frameon=False)

        

    # Adjust where the axes start within the figure. Default value is 10% in from the left and bottom edges. Used to make room for the plot within the figure size (to avoid using bbox_inches='tight' in the savefig-command, as this screws with plot dimensions)
    plt.subplots_adjust(left=options['subplots_adjust'][0], bottom=options['subplots_adjust'][1], right=options['subplots_adjust'][2], top=options['subplots_adjust'][3])


    # If limits for x- and y-axes is passed, sets these.
    if options['xlim'] is not None:
        ax.set_xlim(options['xlim'])

    if options['ylim'] is not None:
        ax.set_ylim(options['ylim'])


    # Add custom text
    if options['text']:
        plt.text(x=options['text'][1][0], y=options['text'][1][1], s=options['text'][0])
    
    return fig, ax



def prettify_thermal_plot_old(fig, ax, options, colour_cycle=('qualitative', 'Dark2_8'), temperatureunit='K', energyunit='eV f.u.$^{-1}$', mode='absolute', scale=1, linestyles=None, labels=None, xpad=None, ypad=None):
    
    # Set sizes of ticks, labes etc.
    ticksize = 30*scale
    labelsize = 30*scale
    legendsize = 30*scale
    titlesize = 30*scale
    
    linewidth = 3*scale
    axeswidth = 3*scale
    majorticklength = 20*scale
    minorticklength = 10*scale

    xpad = 4 if not xpad else xpad
    ypad = 4 if not ypad else ypad
    

    # Set labels on x- and y-axes
    ax.set_xlabel('Temperature [{}]'.format(temperatureunit), size=labelsize, labelpad=xpad)
    if mode == 'absolute':
        ax.set_ylabel('Total energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)
    elif mode == 'difference_plot':
        ax.set_ylabel('Energy difference [{}]'.format(energyunit), size=labelsize, labelpad=ypad)

    elif mode == 'relative':
        ax.set_ylabel('Relative energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)
    
    # Set tick parameters
    ax.tick_params(axis='x', direction='in', which='major', top=True, length=majorticklength, width=axeswidth, pad=xpad)
    ax.tick_params(axis='x', direction='in', which='minor', top=True, length=minorticklength, width=axeswidth)
    
    ax.tick_params(axis='y', direction='in', which='major', right=True, length=majorticklength, width=axeswidth, pad=ypad)
    ax.tick_params(axis='y', direction='in', which='minor', right=True, length=minorticklength, width=axeswidth)
    
    
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    
    if mode == 'absolute':
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
                                   
    elif mode == 'difference_plot':
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))                                   
                                
    
    
    if labels:    
        
        custom_lines = []
        
        if not linestyles:
            linestyles = ['solid', 'dotted', 'dashed']
        
        mod = importlib.import_module("palettable.colorbrewer.%s" % colour_cycle[0])
        palette = getattr(mod, colour_cycle[1]).mpl_colors
        palette = itertools.cycle(palette)
    
        colours = []
        for label in labels:
            colours.append(next(palette))
        
        
        patches = []
        for ind, label in enumerate(labels):
            custom_lines.append(Line2D([0], [0], color=colours[ind], lw=linewidth, ls=linestyles[ind]))
            

        ax.legend(custom_lines, labels, fontsize=labelsize, frameon=False)
    
    
    plt.xticks(fontsize=ticksize, rotation=45)
    plt.yticks(fontsize=ticksize)
    
                                   
    return fig, ax


def prepare_inset_axes(parent_ax, options):
    
    required_options = ['hide_inset_x_labels','hide_inset_x_ticklabels', 'hide_inset_x_ticks', 'rotation_inset_x_ticks', 'hide_inset_y_labels',  'hide_inset_y_ticklabels', 'hide_inset_y_ticks', 'rotation_inset_y_ticks',
        'inset_x_tick_locators', 'inset_y_tick_locators', 'inset_position', 'legend_position']

    default_options = {
        'hide_inset_x_labels': False, # Whether x labels should be hidden
        'hide_inset_x_ticklabels': False,
        'hide_inset_x_ticks': False,
        'rotation_inset_x_ticks': 0,
        'hide_inset_y_labels': False, # whether y labels should be hidden
        'hide_inset_y_ticklabels': False,
        'hide_inset_y_ticks': False,
        'rotation_inset_y_ticks': 0,
        'inset_x_tick_locators': [100, 50], # Major and minor tick locators
        'inset_y_tick_locators': [10, 5],
        'inset_position': [0.1,0.1,0.3,0.3],
        'legend_position': ['upper center', (0.20, 0.90)] # the position of the legend passed as arguments to loc and bbox_to_anchor respectively 
    }
        

    options = update_options(options=options, required_options=required_options, default_options=default_options)


    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    inset_ax = plt.axes([0, 0, 2, 2])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(parent_ax, options['inset_position'])
    inset_ax.set_axes_locator(ip)

    mark_inset(parent_ax, inset_ax, loc1=2, loc2=4, fc='none', ec='black')
    
    inset_ax.xaxis.set_major_locator(MultipleLocator(options['inset_x_tick_locators'][0]))
    inset_ax.xaxis.set_minor_locator(MultipleLocator(options['inset_x_tick_locators'][1]))

    
    inset_ax.yaxis.set_major_locator(MultipleLocator(options['inset_y_tick_locators'][0]))
    inset_ax.yaxis.set_minor_locator(MultipleLocator(options['inset_y_tick_locators'][1]))
   

    
    
    return inset_ax



def update_rc_params(rc_params):
    ''' Update all passed run commands in matplotlib'''
    
    if rc_params:
        for key in rc_params.keys():
            plt.rcParams.update({key: rc_params[key]})

        
        
def update_options(options, required_options, default_options):
    ''' Update all passed options'''
    

    for option in required_options:
        if option not in options.keys():
            options[option] = default_options[option]



    return options



def determine_width(options):
    
    conversion_cm_inch = 0.3937008 # cm to inch
    
    if options['column_type'] == 'single':
        column_width = options['single_column_width']
    elif options['column_type'] == 'double':
        column_width = options['double_column_width']
        
    column_width *= conversion_cm_inch
    
    
    width_ratio = [float(num) for num in options['width_ratio'].split(':')]

    
    width = column_width * width_ratio[0]/width_ratio[1]

    
    return width


def determine_height(options, width):
    
    aspect_ratio = [float(num) for num in options['aspect_ratio'].split(':')]
    
    height = width/(aspect_ratio[0] / aspect_ratio[1])
    
    return height


def scale_figure(options, width, height):
    width = width * options['upscaling_factor'] * options['compress_width']
    height = height * options['upscaling_factor'] * options['compress_height']

    return width, height



def swap_values(dict, key1, key2):

    key1_val = dict[key1]
    dict[key1] = dict[key2]
    dict[key2] = key1_val

    return dict



def generate_colours(palettes):

    # Creates a list of all the colours that is passed in the colour_cycles argument. Then makes cyclic iterables of these. 
    colour_collection = []
    for palette in palettes:
        mod = importlib.import_module("palettable.colorbrewer.%s" % palette[0])
        colour = getattr(mod, palette[1]).mpl_colors
        colour_collection = colour_collection + colour

    colour_cycle = itertools.cycle(colour_collection)


    return colour_cycle