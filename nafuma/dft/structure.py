import math
import re
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import importlib
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from mpl_toolkits.axisartist.axislines import Subplot
from cycler import cycler
import itertools

from ase import Atoms
from ase.io.trajectory import Trajectory
from ase import io
from ase.units import kJ
from ase.eos import EquationOfState
import os
import os.path


def read_eos_data(path, options):
    ''' Reads volume and energy data from a energy-volume run and fits the data to an equation of state. Outputs a list with one pandas DataFrame containing the data points from the DFT-calculations,
    one DataFrame containing the fitted curve data points and one dictionary with equilibrium volume, equilibrium energy and bulk modulus in GPa
    
    path: Path to the folder containing the energ.dat and POSCAR files. energ.dat must have two columns with volumes in the first, energy in the second separated by whitespace. 
    atoms_per_fu: Number of atoms per formula unit. Used to scale the values to be comparable with other calculations that may have a different sized unit cell.
    eos: Type of equation of state to fit to. Same keywords as the ones used in ASE, as it simply calls ASE to fit the equation of state.
    '''

    required_options = ['atoms_per_fu', 'reference', 'eos']

    default_options = {
        'atoms_per_fu': -1, # Scaling factor to output energy per f.u.
        'reference': 0,  # Whether the energy should be relative to some reference energy (typically lowest energy)
        'eos': 'birchmurnaghan', # what type of EoS curve to fit the data to. Options: murnaghan, birch, birchmurnaghan, vinet, pouriertarantola
    }
    

    options = update_options(options=options, required_options=required_options, default_options=default_options)


    # Make paths for the energ.dat and POSCAR files.
    energ_path = os.path.join(path, 'energ.dat')
    poscar_path = os.path.join(path, 'POSCAR')

    # Read POSCAR and calculate the scale factor to give values per formula unit
    at = io.read(poscar_path)
    
    if options['atoms_per_fu'] == -1:
        scale_factor = 1
    else:
        scale_factor = options['atoms_per_fu'] / len(at)
    
    # Get the label
    label = os.path.basename(path)
    
    # Reads the energ.dat file and structures the data into a pandas DataFrame. Then scales the values according to the scale factor.
    dft_df = pd.read_csv(energ_path, delim_whitespace=True, header=None)
    dft_df.columns = ['Configuration', 'Volume', 'Energy']
    dft_df['Energy'] = dft_df['Energy'] * scale_factor
    dft_df['Volume'] = dft_df['Volume'] * scale_factor


    dft_df["Energy"] = dft_df["Energy"] - options['reference'] # subtracts a reference energy if provided. THis value defaults to 0, so will not do anything if not provided.
    
    # Fit data to Equation of State using ASEs EquationOfState object. Makes a DataFrame out of the data points of the fitted curve. Also makes a ditionary of the equilibrium constants, 
    #then packages everything in a list which is returned by the function.
    eos = EquationOfState(dft_df['Volume'].values, dft_df['Energy'].values, eos=options['eos'])
    v0, e0, B = eos.fit()
    
    eos_df = pd.DataFrame(data={'Volume': eos.getplotdata()[4], 'Energy': eos.getplotdata()[5]}) 

    equilibrium_constants = {'v0': v0, 'e0': e0,'B': B/kJ * 1.0e24} 
    
    data = [dft_df, eos_df, equilibrium_constants, label]
    
    return data
    

def read_eos_datas(path, options):


    required_options = ['subset', 'sort_by']

    default_options = {
        'subset': None, # list with directory names of what you want to include
        'sort_by': 'e0', # whether the data should be sorted or not - relevant for bar plots, but also for the order of the entries in the legend in the EoScruve plot
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)
     
    # If a subset of directories is not specified, will create a list of all directories in the path given.
    if not options['subset']:
        dirs = [dir for dir in os.listdir(path) if os.path.isdir(os.path.join(path, dir)) and dir[0] != '.']
    else:
        dirs = options['subset']

    
    datas = []


    # Loop through all subdirectories and reads the data from these. Also appends the name of the directory to the list that is returned from the plot_eos_data() function
    for dir in dirs:
        subdir = os.path.join(path, dir)
        data = read_eos_data(subdir, options)
        datas.append(data)

        
    # Sorts the data if sort is enabled.
    if options['sort_by']:
        datas = sort_data(datas, options['sort_by'])


    return datas
    

def get_summarised_data(path, options):

    datas = read_eos_datas(path=path, options=options)

    summary = []
    for data in datas:
        summary.append([data[3], data[2]['e0'], data[2]['v0'], data[2]['B']])

    df = pd.DataFrame(summary)
    df.columns = ['Label', 'E0', 'V0', 'B']

    emin = df["E0"].min()

    df["dE0"] = df["E0"] - emin

    # Rearranging the columns
    df = df[['Label', 'E0', 'dE0', 'V0', 'B']]

    return df


def plot_eos_data(path, options):
    ''' Plots the data from the energy-volume curve runs. Allows plotting of just the energy-volume curves, a bar plot showing the equilibrium energies or both. 
    
    path: path to where the data is located. It should point to a directory with subdirectories for each structure to be plotted. Inside each of these subdirectories there should be an energ.dat and a POSCAR file.
    atoms_per_fu: Number of atoms per formula unit. Used to scale the values to be comparable with other calculations that may have a different sized unit cell.
    dirs: List of directory names if only a subset of all available datasets is to be plotted. Defaults to None, and will thus get data from all subdirectories.
    eos: Type of equation of state to fit to. Same keywords as the ones used in ASE, as it simply calls ASE to fit the equation of state.
    width: Width of the total figure. Defaults to None, which will again default to width=20.
    height: Height of the total figure. Defaults to None, which will again will default to height= width / phi where phi is the golden ratio.
    dpi: Dots per inch of the figure. Defaults to pyplot's default
    colour_cycles: List of tuples with sets of colours for the palettable colour collection. Defaults to two sets of in total 20 colours. Used for giving different colours to energy-volume curves.
    energyunit: The energy unit. Defaults to eV per formula unit. Only used on the axis labels.
    volumeunit: The volume unit. Defaults to Å^3. Only used on the axis labels.
    xlim: Limits of the x-axes. List of min and max. If mode = both is used, has to contain two lists for each of the plots. As the x-limits for a bar plot is nonesense, should just contain a list with a NoneType. 
    ylim: Limits of the y-axes. List of min and max. If mode = both is used, has to contain two lists for each of the plots.
    sort: Whether or not to sort the data from lowest to highest equilibrium energy. Defaults to True.
    sort_by: What to sort by if sort is enabled. Defaults to e0. Other options: v0 = equilibrium volumes, B = bulk moduli. Alphabetical order sorting is not implemented.
    mode: Determines what to plot. Defaults to energy-volume curves ('curves'). Other options: 'bars', bar-plot of equilibrium energies. 'both', both energy-volume curves and bar plots are plotted side-by-side.
    highlight: Takes a list, either of booleans to highlight certain bars (must be the same length as the number of data sets). Alternatively can contain only names of the datasets to highlight. Defaults to None.'''


    required_options = ['plot_kind', 'highlight', 
                        'reference', 
                        'eos', 'sort_by',
                        'curve_colours',
                        'bar_colours',
                        'marker_cycle',
                        'ylim',
                        'legend_map',
                        'rc_params',
                        'legend']


    default_options = {
        'plot_kind': 'EoScurve', # EoScurve or EoSbars
        'highlight': None, # list with directory names (or Boolean array) of which bars to highlight. Only relevant to EoSbars
        'reference': 0,  # Whether the energy should be relative to some reference energy (typically lowest energy)
        'eos': 'birchmurnaghan', # what type of EoS curve to fit the data to. Options: murnaghan, birch, birchmurnaghan, vinet, pouriertarantola
        'sort_by': 'e0', # whether the data should be sorted or not - relevant for bar plots, but also for the order of the entries in the legend in the EoScruve plot
        'curve_colours': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')], # a set of two colour cycles from the palettable package. Requires many colours for the EoScurve plot
        'bar_colours': [('qualitative', 'Dark2_3'), ('qualitative', 'Pastel2_3')], # set of two colour cycles from the palettable package. Extracts first element of each to make the highlighted and subdued colours respectively. Should be replaced in future by explicit passing of colours
        'marker_cycle': ('o', '*', '^', 'v', 'd', 'H', '8', '>', 'P', 'X'), # marker styles for the EoScurve plot
        'ylim': None, # y-limits (ist)
        'legend': True,
        'legend_map': None, # a dictionary with mappings between the folder names and what should appear in the legend
        'rc_params': None # dictionary of run commands to update plot style
    }


    options = update_options(options=options, required_options=required_options, default_options=default_options)

    # Create path to the data
    datas = read_eos_datas(path=path, options=options)
    
    
    ### PLOT THE ENERGY-VOLUME CURVES
    if options['plot_kind'] == 'EoScurve':
        
        # Fetches a figure and axes object from the prepare_plot() function
        fig, ax = prepare_plot(options=options)
        
        # Make an cyclic iterable of markers to be used for the calculated data points. 
        marker_cycle = itertools.cycle(options['marker_cycle'])
    
        
        # Creates a list of all the colours that is passed in the colour_cycles argument. Then makes cyclic iterables of these. 
        colour_collection = []
        for cycle in options['curve_colours']:
            mod = importlib.import_module("palettable.colorbrewer.%s" % cycle[0])
            colour = getattr(mod, cycle[1]).mpl_colors
            colour_collection = colour_collection + colour

        colour_cycle = itertools.cycle(colour_collection)
        
        labels = []
        colours = []
        markers = []


        # For each of the data sets, extracts the data and plots them.
        for data in datas:
            dft_df, eos_df, label = data[0], data[1], data[3]
            

            # If ylim is passed, only plot those that have a minimum energy below the max ylim parameter
            if options['ylim']:
                plot = True if dft_df["Energy"].min() < options['ylim'][1] else False
            else:
                plot = True

            if plot:    
                if options['label_map']:
                    labels.append(options['label_map'][label])
                
                colours.append(next(colour_cycle))
                markers.append(next(marker_cycle))
        
                dft_df.plot.scatter(x=1, y=2, ax=ax, marker=markers[-1], color=colours[-1], s=20)
                eos_df.plot(x=0, y=1, ax=ax, color=colours[-1], label='_', ls='--')
        
            
        if options['legend']:
            options['legend_content'] = [labels, colours, markers]

        
    ### PLOT THE BAR PLOTS
    elif options['plot_kind'] == 'EoSbars':
        
        # Fetches a figure and axes object from the prepare_plot() function
        fig, ax = prepare_plot(options=options)
        
        e0 = []
        labels = []
        colours = []
        
        # Pick out colour for highlighting (NB! These colours are not passed as arguments, but could be in future)
       
        bar_colours = []
        for cycle in options['bar_colours']:
            mod = importlib.import_module("palettable.colorbrewer.%s" % cycle[0])
            bar_colours.append(getattr(mod, cycle[1]).mpl_colors[0])

       
        # Loops through the datasets, picks out equilibrium volume and labels and sets colours according to the whether the highlight option is used or not.
        for data in datas:

            if options['ylim']:
                plot = True if data[2]['e0'] < options['ylim'][1] else False
            else:
                plot = True

            if plot:

                # Adds 100 if plotting in relative mode. The bases of the bar plots are sunk by 100 during plotting
                adjustment = 100 if options['reference'] != 0 else 100
                print(adjustment)

                e0.append(data[2]['e0']+adjustment)
                print(e0[-1])
                labels.append(options['label_map'][data[3]])
            
                if options['highlight'] is not None:
                    if data[3] in options['highlight']:
                        colours.append(bar_colours[0])
                    else:
                        colours.append(bar_colours[1])

                elif options['highlight'] is not None and type(options['highlight'][0] == str):
                    if labels[-1] in options['highlight']:
                        colours.append(bar_colours[0])
                    else:
                        colours.append(bar_colours[1])

                else:
                    colours.append(bar_colours[0])
                    
        # Makes the bar plot.
        bottom = -100 if options['reference'] != 0 else 0
        plt.bar(range(len(e0)), e0, color=colours, bottom=bottom)
        plt.xticks(range(len(e0)), labels, rotation=90)
        
            
    fig, ax = prettify_plot(fig=fig, ax=ax, options=options)
        
    
    
    
def sort_data(datas, sort_by='e0'):
    ''' Bubble sort algorithm to sort the data sets'''
    
    l = len(datas)
    
    for i in range(0, l):
        for j in range(0, l-i-1):
            if datas[j][2]['{}'.format(sort_by)] > datas[j+1][2]['{}'.format(sort_by)]:
                temp = datas[j]
                datas[j] = datas[j+1]
                datas[j+1] = temp
                
    return datas



def prepare_plot(options={}):
    
    # Reset run commands
    plt.rcdefaults()
    
    # Update run commands if any is passed
    if 'rc_params' in options.keys():
        update_rc_params(options['rc_params'])
        
        
    
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
    'dpi': 600}
    
    options = update_options(options, required_options, default_options)
    
    width = determine_width(options)
    height = determine_height(options, width)
    width, height = scale_figure(options=options, width=width, height=height)
    
    fig, ax = plt.subplots(figsize=(width, height), dpi=options['dpi'])
    
    return fig, ax
    
    
    
        
        
        
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



def prepare_plot_old(width=None, height=None, dpi=None, energyunit='eV', volumeunit=r'Å$^3$', mode='curves', width_ratio=[1, 1], square=True, pad_bottom=None, scale=1, format_params=None):
    '''Prepares pyplot figure and axes objects.'''


    
    linewidth = 3*scale
    axeswidth = 3*scale
    
    plt.rc('lines', linewidth=linewidth)
    plt.rc('axes', linewidth=axeswidth)
    
    
    if square:
        if not width:
            width = 20
            
        height = width
        
        
    else: 
        if not width:
            width = 20
        
        
        if not height:
            golden_ratio = (math.sqrt(5) - 1) / 2
            height = width*golden_ratio
        
    
    
    if mode == 'curves':
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width, height), facecolor='w', dpi=dpi)
    
    
    if mode == 'bars':
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width, height), facecolor='w', dpi=dpi)
        
        
    if mode == 'both':
        
        fig, ax = plt.subplots(1, 2, figsize=(width, height), gridspec_kw={'width_ratios': width_ratio})
        
        
    return fig, ax


def prettify_plot(fig, ax, options):
    '''Prepares pyplot figure and axes objects.'''
    
    required_options = ['plot_kind', 'hide_x_labels', 'hide_y_labels', 'xunit', 'yunit', 'legend_content', 'legend_position', 'x_tick_locators', 'y_tick_locators', 'tick_directions', 'subplots_adjust', 'xlim', 'ylim']

    default_options = {
        'plot_kind': 'EoScurve', # EoScurve or EoSbars
        'hide_x_labels': False, # Whether x labels should be hidden
        'hide_y_labels': False, # whether y labels should be hidden
        'xunit': r'Å$^3$', # The unit of the x-values in the curve plot
        'yunit': r'eV f.u.$^{-1}$', # The unit of the y-values in the curve and bar plots
        'xlim': None,
        'ylim': None,
        'legend_content': None,
        'legend_position': ['upper center', (1.10, 0.90)], # the position of the legend passed as arguments to loc and bbox_to_anchor respectively
        'x_tick_locators': [10, 5], # Major and minor tick locators
        'y_tick_locators': [.1, .05], # Major and minor tick locators
        'tick_directions': 'in', # in or out
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9]
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    
    if options['plot_kind'] == 'EoScurve':
   
        # Set labels on x- and y-axes
        ax.set_xlabel('Volume [{}]'.format(options['xunit']))
        
        if not options['hide_y_labels']:
            ax.set_ylabel('Energy [{}]'.format(options['yunit']))
        else:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)


        ax.xaxis.set_major_locator(MultipleLocator(options['x_tick_locators'][0]))
        ax.xaxis.set_minor_locator(MultipleLocator(options['x_tick_locators'][1]))

        ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
        ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))
        
        if ax.get_legend():
            ax.get_legend().remove()
        

        if options['legend']:
            labels = options['legend_content'][0]
            colours = options['legend_content'][1]
            markers = options['legend_content'][2]
            
            entries = []
            
            for i in range(len(options['legend_content'][0])):
                entries.append(mlines.Line2D([], [], label=labels[i], color=colours[i], marker=markers[i], linestyle='None'))
                
            
            fig.legend(handles=entries, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1], frameon=False)
            

  
        
    if options['plot_kind'] == 'EoSbars':
        
        if not options['hide_y_labels']:
            ax.set_ylabel('Energy [{}]'.format(options['yunit']))
        
        ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
        ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))

        ax.tick_params(axis='x', which='minor', bottom=False, top=False)
        



    # Adjust where the axes start within the figure. Default value is 10% in from the left and bottom edges. Used to make room for the plot within the figure size (to avoid using bbox_inches='tight' in the savefig-command, as this screws with plot dimensions)
    plt.subplots_adjust(left=options['subplots_adjust'][0], bottom=options['subplots_adjust'][1], right=options['subplots_adjust'][2], top=options['subplots_adjust'][3])


    # If limits for x- and y-axes is passed, sets these.
    if options['xlim'] is not None:
        ax.set_xlim(options['xlim'])

    if options['ylim'] is not None:
        ax.set_ylim(options['ylim'])
        
    
    return fig, ax



def prettify_plot_old(fig, ax, energyunit='eV', volumeunit=r'Å$^3$', mode='curves', legend_content=None, pad_bottom=None, scale=1, hide_ylabels=False, xpad=None, ypad=None):
    '''Prepares pyplot figure and axes objects.'''
    
    # Set sizes of ticks, labes etc.
    ticksize = 30*scale
    labelsize = 30*scale
    legendsize = 15*scale
    titlesize = 30*scale
    
    linewidth = 3*scale
    axeswidth = 3*scale
    markersize = 15*scale
    majorticklength = 20*scale
    minorticklength = 10*scale

    xpad = 4 if not xpad else xpad
    ypad = 4 if not ypad else ypad

    
    if mode == 'curves':
   
        # Set labels on x- and y-axes
        ax.set_xlabel('Volume [{}]'.format(volumeunit), size=labelsize, labelpad=xpad)
        
        if not hide_ylabels:
            ax.set_ylabel('Energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)
        else:
            ax.set_ylabel('')

        # Set tick parameters
        ax.tick_params(axis='both', direction='in', which='major', length=majorticklength, width=axeswidth, right=True, top=True, labelsize=ticksize)
        ax.tick_params(axis='both', direction='in', which='minor', length=minorticklength, width=axeswidth, right=True, top=True, labelsize=ticksize)

        ax.tick_params(axis='x', pad=xpad)
        ax.tick_params(axis='y', pad=ypad)

        if hide_ylabels:
            ax.tick_params(labelleft=False)

        plt.xticks(fontsize=ticksize)
        plt.yticks(fontsize=ticksize)
        

        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))

        ax.yaxis.set_major_locator(MultipleLocator(.1))
        ax.yaxis.set_minor_locator(MultipleLocator(.05))
        

        ax.get_legend().remove()
        if legend_content:
            patches = []
            labels = legend_content[0]
            colours = legend_content[1]
            markers = legend_content[2]
            
            entries = []
            
            for ind, label in enumerate(legend_content[0]):
                entries.append(mlines.Line2D([], [], color=colours[ind], marker=markers[ind], linestyle='None',
                          markersize=markersize, label=labels[ind]))
                
                #patches.append(mpatches.Patch(color=colours[ind], label=labels[ind]))
        
            
            fig.legend(handles=entries, loc='upper center', bbox_to_anchor=(1.10, 0.90), fontsize=legendsize, frameon=False)
            
        if pad_bottom is not None:
            bigax = fig.add_subplot(111) 
            bigax.set_facecolor([1,1,1,0])
            bigax.spines['top'].set_visible(False)
            bigax.spines['bottom'].set_visible(True)
            bigax.spines['left'].set_visible(False)
            bigax.spines['right'].set_visible(False)
            bigax.tick_params(labelcolor='w', color='w', direction='in', top=False, bottom=True, left=False, right=False, labelleft=False, pad=pad_bottom)
        
    if mode == 'bars':
        
        
        ax.tick_params(axis='both', direction='in', which='major', length=majorticklength, width=axeswidth, right=True, top=True)
        ax.tick_params(axis='both', direction='in', which='minor', length=minorticklength, width=axeswidth, right=True, top=True)
        
        if not hide_ylabels:
            ax.set_ylabel('Energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)
        
        ax.yaxis.set_major_locator(MultipleLocator(.1))
        ax.yaxis.set_minor_locator(MultipleLocator(.05))

        ax.tick_params(axis='x', pad=xpad)
        ax.tick_params(axis='y', pad=ypad)
        
        plt.xticks(fontsize=ticksize)
        plt.yticks(fontsize=ticksize)
        
        if pad_bottom is not None:
            bigax = fig.add_subplot(111) 
            bigax.set_facecolor([1,1,1,0])
            bigax.spines['top'].set_visible(False)
            bigax.spines['bottom'].set_visible(True)
            bigax.spines['left'].set_visible(False)
            bigax.spines['right'].set_visible(False)
            bigax.tick_params(labelcolor='w', color='w', direction='in', top=False, bottom=True, left=False, right=False, labelleft=False, pad=pad_bottom)
        
    if mode == 'both':
        
        # Set labels on x- and y-axes
        ax[0].set_xlabel('Volume [{}]'.format(volumeunit), size=labelsize, labelpad=xpad)
        ax[0].set_ylabel('Energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)

        # Set tick parameters
        ax[0].tick_params(axis='both', direction='in', which='major', length=majorticklength, width=axeswidth, right=True, left=True, top=True, labelsize=ticksize)
        ax[0].tick_params(axis='both', direction='in', which='minor', length=minorticklength, width=axeswidth, right=True, left=True, top=True, labelsize=ticksize)

        ax[0].tick_params(axis='x', pad=xpad)
        ax[0].tick_params(axis='y', pad=ypad)
       
        ax[0].xaxis.set_major_locator(MultipleLocator(10))
        ax[0].xaxis.set_minor_locator(MultipleLocator(5))

        ax[0].yaxis.set_major_locator(MultipleLocator(.1))
        ax[0].yaxis.set_minor_locator(MultipleLocator(.05))
        
        plt.xticks(fontsize=ticksize)
        plt.yticks(fontsize=ticksize)


        ax[1].yaxis.set_major_locator(MultipleLocator(.2))
        ax[1].yaxis.set_minor_locator(MultipleLocator(.1))
        ax[1].yaxis.set_label_position('right')
        ax[1].yaxis.tick_right()
        ax[1].set_ylabel('Energy [{}]'.format(energyunit), size=labelsize, ypad=ypad)
        ax[1].tick_params(axis='both', direction='in', which='major', length=majorticklength, width=axeswidth, left=True, right=True, top=True)
        ax[1].tick_params(axis='both', direction='in', which='minor', length=minorticklength, width=axeswidth, left=True, right=True, top=True)

        ax[1].tick_params(axis='x', pad=xpad)
        ax[1].tick_params(axis='y', pad=ypad)
        
        
        plt.xticks(fontsize=ticksize)
        plt.yticks(fontsize=ticksize)
        
    return fig, ax



def parabola(V, a, b, c):
    """parabola polynomial function

    this function is used to fit the data to get good guesses for
    the equation of state fits

    a 4th order polynomial fit to get good guesses for
    was not a good idea because for noisy data the fit is too wiggly
    2nd order seems to be sufficient, and guarantees a single minimum"""


    E = (a * V**2) + (b * V) + c

    return E


def murnaghan(V, E0, V0, B0, BP):
    'From PRB 28,5480 (1983'

    E = E0 + ((B0 * V) / BP) * (((V0 / V)**BP) / (BP - 1) + 1) - ((V0 * B0) / (BP - 1))
    return E


def birch(V, E0, V0, B0, BP):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    paper downloaded from Web

    case where n=0
    """

    E = (E0 +
         9 / 8 * B0 * V0 * ((V0 / V)**(2 / 3) - 1)**2 +
         9 / 16 * B0 * V0 * (BP - 4) * ((V0 / V)**(2 / 3) - 1)**3)
    return E


def birchmurnaghan(V, E0, V0, B0, BP):
    """
    BirchMurnaghan equation from PRB 70, 224107
    Eq. (3) in the paper. Note that there's a typo in the paper and it uses
    inversed expression for eta.
    """

    eta = (V0 / V)**(1 / 3)

    E = E0 + 9 * B0 * V0 / 16 * (eta**2 - 1)**2 * (6 + BP * (eta**2 - 1) - 4 * eta**2)

    return E


def vinet(V, E0, V0, B0, BP):
    'Vinet equation from PRB 70, 224107'

    eta = (V / V0)**(1 / 3)

    E = (E0 + 2 * B0 * V0 / (BP - 1)**2 *
         (2 - (5 + 3 * BP * (eta - 1) - 3 * eta) *
          np.exp(-3 * (BP - 1) * (eta - 1) / 2)))
    
    return E

def pouriertarantola(V, E0, V0, B0, BP):
    'Pourier-Tarantola equation from PRB 70, 224107'

    eta = (V / V0)**(1 / 3)
    squiggle = -3 * np.log(eta)

    E = E0 + B0 * V0 * squiggle**2 / 6 * (3 + squiggle * (BP - 2))
    return E



def get_initial_guesses(volume, energy):

    p = np.polyfit(volume, energy, deg=2)

    a, b, c = p[0], p[1], p[2]

    # Estimated from dE/dV = 2aV0 + b => V0 = -b / 2a
    v0 = -b / (2*a)

    # Estimated by evaluating a parabola with a, b and c values at V = V0
    e0 = parabola(v0, a, b, c)

    # Estimated form B0 ~ V0 * d^2E / dV^2. d^2E / dV^2 = 2a. 
    b0 = 2 * a * v0

    # Just a reasonable starting value
    bp = 4


    return [e0, v0, b0, bp]



def fit_eos_curve(volume, energy, p0, eos):

    eos_dict = {'murnaghan': murnaghan, 'birch': birch, 'birchmurnaghan': birchmurnaghan, 'vinet': vinet, 'pouriertarantola': pouriertarantola}

    func = eos_dict[eos]

    popt, pcov = curve_fit(func, volume, energy, p0)

    E0, V0, B0, BP  = popt[0], popt[1], popt[2], popt[3]

    return [E0, V0, B0, BP]




def get_plotdata(volume, energy, equilibrium_values, eos):
 
    eos_dict = {'murnaghan': murnaghan, 'birch': birch, 'birchmurnaghan': birchmurnaghan, 'vinet': vinet, 'pouriertarantola': pouriertarantola}

    V = np.linspace(volume.min(), volume.max(), 100)

    E0, V0, B0, BP = equilibrium_values[0], equilibrium_values[1], equilibrium_values[2], equilibrium_values[3]

    print(E0, V0, B0, BP)

    func = eos_dict[eos]

    print(func)

    E = func(V, E0, V0, B0, BP)

    return E, V


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



def get_equilibrium_data(path, atoms_per_formula_unit, eos=None):


    if not eos:
        eos = 'murnaghan'


    dirs = [os.path.join(path, dir) for dir in os.listdir(path)]



    data = []

    for dir in dirs:
        atoms, atom_num, atoms_dict = get_atoms(os.path.join(dir, 'POSCAR'))
        scaling_factor = sum(atom_num) / atoms_per_formula_unit

        label = dir.split('/')[-1]

        dft_df = pd.read_csv(os.path.join(dir, 'energ.dat'), header=None, delim_whitespace=True)
        dft_df.columns = ['Volume', 'Energy']

        volume = dft_df["Volume"].to_numpy() / scaling_factor
        energy = dft_df["Energy"].to_numpy() / scaling_factor

        p0 = get_initial_guesses(volume, energy)

        equilibrium_constants = fit_eos_curve(volume, energy, p0, eos)
        e0, v0, b0, bp = equilibrium_constants[0], equilibrium_constants[1], equilibrium_constants[2], equilibrium_constants[3]

        data.append([label, e0, v0, b0/kJ*1e24, bp])


    df = pd.DataFrame(data)
    df.columns = ['Label', 'E0', 'V0', 'B0', 'Bp']
    df.sort_values(by='E0', ascending=True, inplace=True)
    df.reset_index(inplace=True)

    E_min = df['E0'].min()

    df['dE'] = df['E0'] - E_min

    df = df[['Label', 'E0', 'dE', 'V0', 'B0', 'Bp']]


    return df




