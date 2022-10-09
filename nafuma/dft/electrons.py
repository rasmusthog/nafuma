import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import linecache
import nafuma.dft as dft
import nafuma.auxillary as aux
import nafuma.plotting as btp

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import importlib
import mpl_toolkits.axisartist as axisartist
from cycler import cycler
import itertools
import matplotlib.patches as mpatches




def count_electrons(pdos, orbital, interval=None, r=None, scale=None):
    ''' Counts electrons the specified oribtals from a projected density of states DataFrame. Interval can be specified, as well as a scaling factor and whether the number should be rounded. 
    Inputs:
    dos: either an individual DOS (as read from read_pdos()), or a list of DOSes. If a single DataFrame is passed, it will be appended to a list
    orbital: list of which orbitals to count the electrons from
    interval: a list specifying where to start counting from (lower limit) to where to stop counting (upper limit) in eV
    r: Number of decimals points the number should be rounded to
    scale: A scaling factor to scale the number of electrons to a desired size, e.g. if you have a set containing two atoms per unit cell and you want to know how many electrons per atom there is

    Output:
    nelec: The total number of electrons given your choices
    nelec_dos: A list where each element is the total number of electrons per DOS passed (e.g. you pass three PDOS from three individual atoms, then you will get total electron count per atom)
    nelec_orbitals: A list of lists, where each list contains the number of electrons per orbital specified (e.g. you pass three PDOS from three individual atoms, you will get three lists each containing electron count per orbital specified)
    '''

    

    if not type(orbital) == list:
        orbital = [orbital]

    if not type(pdos) == list:
        pdos = [pdos]

    
    nelec = 0
    nelec_per_dos = []
    nelec_per_orbital = []
    
    for d in pdos:

        energy = d["Energy"]

        nelec_orbitals = []

        for o in orbital:
            orbital_dos = d[o]
            dE = (energy.max()-energy.min()) / len(energy)


            if not interval:
                interval = [energy.min(), energy.max()]

            emin, emax = interval[0], interval[1]


            nelec_orbital = 0
            for od, e in zip(orbital_dos, energy):
                if e > emin and e < emax:
                    nelec_orbital += np.abs(od)*dE
                    #print(nelec_orbital)

            nelec += nelec_orbital
            nelec_orbitals.append(nelec_orbital)

        
        # Scale the values if specified
        if scale:
            
            for ind, nelec_orbital in enumerate(nelec_orbitals):
                nelec_orbitals[ind] = nelec_orbital * scale


        # If rounding is specified, does so to the electron count per DOS and the electron count per orbital
        if r:
            # First sums the electron count per orbital, and then round this number
            nelec_dos = np.round(sum(nelec_orbitals), r)

            # Then each individual orbital electron count
            for ind, nelec_orbital in enumerate(nelec_orbitals):
                nelec_orbitals[ind] = np.round(nelec_orbital, r)

        # If no rounding is specified, just adds the electron count per orbital together
        else:
            nelec_dos = sum(nelec_orbitals)

        # Appends the total electron count for this DOS to the list of all individual DOS electron count
        nelec_per_dos.append(nelec_dos)

        # Appends the list of orbital electron counts to the list of all the individual DOS orbital electron count (phew...)
        nelec_per_orbital.append(nelec_orbitals)


    # The total electron count is then scaled in the end. At this point the other values will have been scaled already
    if scale:
        nelec = nelec * scale

    # And lastly round if this is specified. Again, the electron counts in the lists are already rounded so they don't have to be rounded again 
    if r:
        nelec = np.round(nelec, r)

    return nelec, [nelec_per_dos, nelec_per_orbital]



def integrate_coop(coopcar, interval=None, r=None, scale=None, interactions=None, kind='individual', up=True, down=True, collapse=False):
    ''' As of now does not support not passing in interactions. Very much copy and paste from the plotting function - not every choice here might make sense for integration of COOP'''

    coopcar, coop_interactions = dft.io.read_coop(coopcar, collapse=collapse)

    # If interactions has been specified
    if interactions:

        # Make interactions into a list of lists for correct looping below
        if type(interactions[0]) != list:
            interactions_list = [interactions]
        else:
            interactions_list = interactions

        for ind, interactions in enumerate(interactions_list):

            # Determine which columns to integrate if collapse is enabled
            if collapse:
                to_integrate = [2*(i-1)+3 for i in interactions]


                # Make mean column for integration if mean mode is enabeld (is this sensible to include?)
                if kind == 'avg' or kind == 'average' or kind == 'mean':
                    coopcar["mean"] = coopcar.iloc[:, to_integrate].mean(axis=1)
                    to_integrate = [coopcar.columns.get_loc('mean')]
            
            # Determine which columns to integrate if collapse is disabled and both up and down should be plotted
            elif up and down:
                to_integrate_up = [2*(i-1)+3 for i in interactions]
                to_integrate_down = [2*(i-1)+5 +2*len(coop_interactions) for i in interactions]
                to_integrate = to_integrate_up + to_integrate_down

                if kind == 'avg' or kind == 'average' or kind == 'mean':
                    coopcar["mean_up"] = coopcar.iloc[:, to_integrate_up].mean(axis=1)
                    coopcar["mean_down"] = coopcar.iloc[:, to_integrate_down].mean(axis=1)
                    to_integrate = [coopcar.columns.get_loc('mean_up'), coopcar.columns.get_loc('mean_down')]

            # Determine which columns to plot if collapse is disabled and only up should be plotted
            elif up:
                to_integrate = [2*(i-1)+3 for i in interactions]

                if kind == 'avg' or kind == 'average' or kind == 'mean':
                    coopcar["mean_up"] = coopcar.iloc[:, to_integrate].mean(axis=1)
                    to_integrate = [coopcar.columns.get_loc('mean_up')]


             # Determine which columns to plot if collapse is disabled and only down should be plotted
            elif down:
                to_integrate = [2*(i-1)+5 +2*len(coop_interactions) for i in interactions]

                if kind == 'avg' or kind == 'average' or kind == 'mean':
                    coopcar["mean_down"] = coopcar.iloc[:, to_integrate].mean(axis=1)
                    to_integrate = [coopcar.columns.get_loc('mean_down')]



        bonding = 0
        antibonding = 0
        bonding_interactions = []
        antibonding_interactions = []
        difference_interactions = []
        percentage_bonding_interactions = []


        for integrate in to_integrate:

            bonding_interaction = 0
            antibonding_interaction = 0

            coop = coopcar.iloc[:, integrate]

            energy = coopcar["Energy"]
            dE = (energy.max()-energy.min()) / len(energy)

            # Sets interval to everything below the Fermi-level by default if not specified in function call
            if not interval:
                interval = [energy.min(), 0]

            emin, emax = interval[0], interval[1]


            for c, e in zip(coop, energy):
                if e > emin and e < emax:
                    if c > 0:
                        bonding_interaction += c*dE
                    elif c < 0:
                        antibonding_interaction += np.abs(c)*dE

            
            bonding += bonding_interaction
            antibonding += antibonding_interaction

            difference_interaction = bonding_interaction - antibonding_interaction
            percentage_bonding_interaction = bonding_interaction / (bonding_interaction + antibonding_interaction) * 100

            if scale:
                bonding_interaction = bonding_interaction * scale
                antibonding_interaction = antibonding_interaction * scale
                difference_interaction = difference_interaction * scale

            if r:
                bonding_interaction = np.round(bonding_interaction, r)
                antibonding_interaction = np.round(antibonding_interaction, r)
                difference_interaction = np.round(difference_interaction, r)
                percentage_bonding_interaction = np.round(percentage_bonding_interaction, r)

            bonding_interactions.append(bonding_interaction)
            antibonding_interactions.append(antibonding_interaction)
            difference_interactions.append(difference_interaction)
            percentage_bonding_interactions.append(percentage_bonding_interaction)

        difference = bonding - antibonding
        percentage_bonding = (bonding/(bonding+antibonding)) * 100

        if scale:
            bonding = bonding * scale
            antibonding = antibonding * scale
            difference = difference * scale

        if r:
            bonding = np.round(bonding, r)
            antibonding = np.round(antibonding, r)
            difference = np.round(difference, r)
            percentage_bonding = np.round(percentage_bonding, r)

    return [bonding, antibonding, difference, percentage_bonding], [bonding_interactions, antibonding_interactions, difference_interactions, percentage_bonding_interactions]




def plot_pdos(data: dict, options={}):

    default_options = {
        'xlabel': 'Energy',                     'xunit': 'eV',      'xlim': None,   'x_tick_locators': None,
        'ylabel': 'Partial density of states',  'yunit': 'arb.u.',  'ylim': None,   'y_tick_locators': None,
        'mark_fermi_level': True,           # Adds a dashed line to mark the Fermi-level
        'flip_axes':        False,          # Flips x- and y-axes
        'plot_indices':     [],             # List which indices to plot. If options["sum_atoms"] == True, this needs to be a list of lists, each specifying the index of a given atom
        'plot_atoms':       [],             # List of which atoms to plot. Only used if options["sum_atoms"] == True.
        'plot_orbitals':    [],             # List of which orbitals to plot. If options["sum_atoms"] == True, this needs to be a list of lists, each specifying the orbitals of a given atom
        'atom_colours':     [],             # Colours of each atom. Should be a colour for each atom, only in use if options["sum_atoms"] == True.
        'orbital_colours':  [],             # Colours of each orbital. The list should always correspond to the shape of options["plot_orbitals"].  
        'fill':             False,
        'fig':              None,           # Matplotlib Figure object
        'ax':               None,           # Matplotlib Axes object
    }

    options = aux.update_options(options=options, default_options=default_options)

    if 'axes_flipped' not in options.keys():
        options['axes_flipped'] = False

    data = dft.io.read_pdos(data=data, options=options)


   
    if not options['fig'] and not options['ax']:
        fig, ax = btp.prepare_plot(options=options)
    else:
        fig, ax = options['fig'], options['ax']


    # If options['sum_atoms'] == True
    if isinstance(data['pdos'], dict):

        # Populate the plot_atoms and plot_orbitals lists if they are not passed. Defaults to showing everything
        if not options['plot_atoms']:
            options['plot_atoms'] = data['atoms']['specie']

        if not options['plot_orbitals']:
            for atom in options['plot_atoms']:
                options['plot_orbitals'].append([])

        # This is to fill in each orbital list for each atom. This is in case options['plot_orbitals'] is passes, but one or more of the atoms lack colours
        for i, atom in enumerate(options['plot_atoms']):
            if not options['plot_orbitals'] or not options['plot_orbitals'][i]:
                options['plot_orbitals'][i] = [orbital for orbital in data['pdos'][atom].columns if 'Energy' not in orbital]

        # Populate the atom_colours and orbital_colours. Defaults to same colour for all orbitals of one specie.
        if not options['atom_colours']:
            options['palettes'] = [('qualitative', 'Dark2_8')]
            colour_cycle = generate_colours(options=options)

            for atom in options['plot_atoms']:
                options['atom_colours'].append(next(colour_cycle))
        
        if not options['orbital_colours']:
            for i, atom in enumerate(options['plot_orbitals']):
                options['orbital_colours'].append([]) # Make list for specific atom
                for orbital in atom:
                    options['orbital_colours'][i].append(options['atom_colours'][i])
             

        for i, atom in enumerate(options['plot_atoms']):
                
            if not options['plot_orbitals'] or not options['plot_orbitals'][i]:
                options['plot_orbitals'][i] = [orbital for orbital in data['pdos'][atom].columns if 'Energy' not in orbital]

            x = 'Energy'
            y = options['plot_orbitals'][i]

            if options['flip_axes']:
                for j, orbital in enumerate(options['plot_orbitals'][i]):

                    if options['fill']:
                        ax.fill_betweenx(y=data['pdos'][atom]['Energy'], x1=data['pdos'][atom][orbital], x2=0, color=options['orbital_colours'][i][j], alpha=0.5, ec=(0,0,0,0))
                    else:
                        ax.plot(data['pdos'][atom][orbital], data['pdos'][atom]['Energy'], color=options['orbital_colours'][i][j])



            else:
                data['pdos'][atom].plot(x=x, y=y, color=options['orbital_colours'][i], ax=ax)

        #print(options['plot_orbitals'], options['orbital_colours'])


    if options['flip_axes']:
        
        if not options['axes_flipped']:
            options = aux.swap_values(options=options, 
                                                key1=['xlim', 'xunit', 'xlabel', 'x_tick_locators'], 
                                                key2=['ylim', 'yunit', 'ylabel', 'y_tick_locators']
                                            )

            options['axes_flipped'] = True # 

        ax.axvline(x=0, c='black')

        if options['mark_fermi_level']:
            ax.axhline(y=0, c='black', ls='--')

    
    else:
        ax.axhline(y=0, c='black')
        ax.axvline(x=0, c='black', ls='--')

    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)




    #elif isinstance(data['pdos'], list):
    #       if not options['plot_atoms']:
    #       options['plot_atoms'] = data['atoms']['specie']
    #   
    #   if not options['plot_indices']:
    #       for plot_specie in options['plot_atoms']:
    #           for i, doscar_specie in enumerate(data['atoms']['specie']):
    #               if plot_specie == doscar_specie:
    #                   options['plot_indices'].append([k for k in range(data['atoms']['number'][i])])


    return None
    




def plot_partial_dos_legacy(data: dict, options={}):


    required_options = ['atoms', 'orbitals', 'up', 'down', 'sum_atoms', 'collapse_spin', 'sum_orbitals', 'palettes', 'colours', 'fig', 'ax']

    default_options = {
        'atoms': None,
        'orbitals': None,
        'up': True,
        'down': True,
        'sum_atoms': True,
        'collapse_spin': False,
        'sum_orbitals': False,
        'palettes': [('qualitative', 'Dark2_8')],
        'colours': None,
        'fig': None,
        'ax': None


    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)

    if not options['ax'] and not options['fig']:
        fig, ax = btp.prepare_plot(options)
    else:
        fig, ax = options['fig'], options['ax']
    
    species, *_ = dft.io.get_atoms(data['poscar'])

    pdos, options['dos_info'] = dft.io.read_pdos(data=data, options=options) # Extract projected DOS from DOSCAR, decomposed on individual atoms and orbitals Should yield list of N DataFrames where N is number of atoms in POSCAR


    if not options['orbitals']:
            options['orbitals'] = ['s', 'p1', 'p2', 'p3', 'd1', 'd2', 'd3', 'd4', 'd5'] if not options['sum_orbitals'] else ['s', 'p', 'd'] 

    if not options['colours']:
        colour_cycle = generate_colours(options=options)
# 
        #colours = []
        #for orbital in options['orbitals']:
            # colours.append(next(colour_cycle))
# 
    # else:
        # colours = options['colours']

    elif not isinstance(options['colours'], list):
        new_colours = []
        for atom in options['atoms']:
            new_colours.append([options['colours']])

        options['colours'] = new_colours


    print(options['colours'])
    

    if not isinstance(options['orbitals'][0], list):
        new_orbitals = []
        for atom in options['atoms']:
            new_orbitals.append([options['orbitals']])

        options['orbitals'] = new_orbitals


    if options['atoms']:
        for i, atom in enumerate(options['atoms']):

            if options['sum_atoms']:
                for ind, specie in enumerate(species):
                    if specie == atom:
                        atom_index = ind
            else:
                atom_index = atom-1

            
            for j, orbital in enumerate(options['orbitals'][i]):
                
                colour = options['colours'][i][j]

                if options['dos_info']['spin_polarised']:
                    if options['up']:
                        pdos[atom_index].plot(x='Energy', y=orbital+'_u', ax=ax, c=colour)

                    if options['down']:
                        pdos[atom_index].plot(x='Energy', y=orbital+'_d', ax=ax, c=colour)
                else:
                    pdos[atom_index].plot(x='Energy', y=orbital, ax=ax, c=colour)



    btp.adjust_plot(fig=fig, ax=ax, options=options)

    return [pdos, ax, fig]



def get_pdos_indices(poscar, atoms):

    species, atom_num, atoms_dict = dft.io.get_atoms(poscar)




def get_pdos(doscar='DOSCAR', nedos=301, headerlines=6, spin=True, adjust=True, manual_adjust=None):
    
    lines = dft.io.open_doscar(doscar)
    
    number_of_atoms = dft.io.get_number_of_atoms(doscar)

    if adjust:
        e_fermi = dft.io.get_fermi_level(doscar) if not manual_adjust else manual_adjust
    else:
        e_fermi = 0
    
    pdos = []
    
    columns_non_spin = ["Energy", "s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z2-r2", "d_xz", "d_x2-y2"]
    columns_spin = ["Energy", "s_up", "s_down", "p_y_up", "p_y_down", "p_z_up", "p_z_down", "p_x_up", "p_x_down", "d_xy_up", "d_xy_down", "d_yz_up", "d_yz_down", 
                        "d_z2-r2_up", "d_z2-r2_down", "d_xz_up", "d_xz_down", "d_x2-y2_up", "d_x2-y2_down"]

    up = ['s_up', "p_y_up",  "p_z_up",  "p_x_up", "d_xy_up", "d_yz_up", "d_z2-r2_up", "d_xz_up", "d_x2-y2_up"]
    down = ['s_down', "p_y_down", "p_z_down", "p_x_down", "d_xy_down", "d_yz_down", "d_z2-r2_down", "d_xz_down", "d_x2-y2_down"]
    total = ["s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z2-r2", "d_xz", "d_x2-y2"]
    
    for i in range(1,number_of_atoms+1):
        atom_dos = []
        
        for j in range(headerlines+(nedos*i)+i,nedos+headerlines+(nedos*i)+i):
            line = lines[j].strip()
            values = line.split()
    
            for ind, value in enumerate(values):
                values[ind] = float(value)
        
            values[0] = values[0] - e_fermi
            atom_dos.append(values)
        
        
        atom_df = pd.DataFrame(data=atom_dos, columns=columns_non_spin) if spin==False else pd.DataFrame(data=atom_dos, columns=columns_spin)
        
        if spin==True:
            atom_df[["s_down"]] = -atom_df[["s_down"]]
            atom_df[["p_y_down"]] = -atom_df[["p_y_down"]]
            atom_df[["p_z_down"]] = -atom_df[["p_z_down"]]
            atom_df[["p_x_down"]] = -atom_df[["p_x_down"]]
            atom_df[["d_xy_down"]] = -atom_df[["d_xy_down"]]
            atom_df[["d_yz_down"]] = -atom_df[["d_yz_down"]]
            atom_df[["d_z2-r2_down"]] = -atom_df[["d_z2-r2_down"]]
            atom_df[["d_xz_down"]] = -atom_df[["d_xz_down"]]
            atom_df[["d_x2-y2_down"]] = -atom_df[["d_x2-y2_down"]]

            atom_df = atom_df.assign(total_up = atom_df[up].sum(axis=1))
            atom_df = atom_df.assign(total_down = atom_df[down].sum(axis=1))

        elif spin==False:
            atom_df = atom_df.assign(total = atom_df[total].sum(axis=1))
        
        pdos.append(atom_df)
        
    return pdos



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


def prepare_dos_plot(width=None, height=None, square=True, dpi=None, colour_cycle=('qualitative', 'Dark2_8'), energyunit='eV', dosunit='arb. u.', scale=1, pdos=None):
    
    if not pdos:
        linewidth = 3*scale
    else:
        linewidth = 3

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

def prettify_dos_plot(fig, ax, options):
    

    required_options = ['plot_kind', 'flip_xy', 'hide_x_labels', 'hide_y_labels',  'xlabel', 'ylabel', 'xunit', 'yunit', 'xlim', 'ylim', 'x_tick_locators', 'y_tick_locators', 'y_tick_format', 'x_tick_format', 'hide_x_ticks', 'hide_y_ticks', 'hide_x_ticklabels', 'hide_y_ticklabels',
                        'colours', 'palettes',  'title', 'legend', 'labels', 'label_colours', 'legend_position', 'legend_ncol', 'subplots_adjust', 'text']

    default_options = {
        'plot_kind': 'PDOS', # DOS/PDOS/COOP/COHP
        'flip_xy': False,
        'hide_x_labels': False, # Whether x labels should be hidden
        'hide_x_ticklabels': False,
        'hide_x_ticks': False,
        'hide_y_labels': False, # whether y labels should be hidden
        'hide_y_ticklabels': False,
        'hide_y_ticks': False,
        'xlabel': 'Energy',
        'ylabel': 'DOS',
        'xunit': r'eV', # The unit of the x-values in the curve plot
        'yunit': r'a.u.', # The unit of the y-values in the curve and bar plots
        'xlim': None,
        'ylim': None,
        'x_tick_locators': [1, .5], # Major and minor tick locators
        'y_tick_locators': [1, .5],
        'y_tick_format': None,
        'x_tick_format': None,
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'title': None,
        'legend': True,
        'labels': None,
        'label_colours': None,
        'legend_position': ['upper center', (0.20, 0.90)], # the position of the legend passed as arguments to loc and bbox_to_anchor respectively 
        'legend_ncol': 1,
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9],
        'text': None
    }


    if 'plot_kind' in options.keys():
        if 'ylabel' not in options.keys():
            if options['plot_kind'] == 'DOS':
                options['ylabel'] = 'DOS'
            elif options['plot_kind'] == 'PDOS':
                options['ylabel'] = 'PDOS'
            elif options['plot_kind'] == 'COOP':
                options['ylabel'] = 'COOP'
            elif options['plot_kind'] == 'COHP':
                options['ylabel'] = 'COHP'
            


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

    # Change format of axis tick labels if specified:


           
    # Set title
    if options['title']:
        ax.set_title(options['title'])


    if options['y_tick_format']:
        ax.yaxis.set_major_formatter(FormatStrFormatter(options['y_tick_format']))
    if options['x_tick_format']:
        ax.xaxis.set_major_formatter(FormatStrFormatter(options['x_tick_format']))


    # Create legend

    if ax.get_legend():
                ax.get_legend().remove()

            
    if options['legend'] and options['labels']:


        # Generate colours
        if not options['colours']:
            colour_cycle = generate_colours(palettes=options['palettes'])
            
            colours = []
            for label in options['labels']:
                colours.append(next(colour_cycle))
        
        
        else:
            colours = options['colours']

        if options['label_colours']:
            colours = options['label_colours']

        # Create legend
        patches = []
        for i, label in enumerate(options['labels']):
            patches.append(mpatches.Patch(color=colours[i], label=label))

        print(options['legend_ncol'])
    
        ax.legend(handles=patches, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1], frameon=False, ncol=options['legend_ncol'])

        

    # Adjust where the axes start within the figure. Default value is 10% in from the left and bottom edges. Used to make room for the plot within the figure size (to avoid using bbox_inches='tight' in the savefig-command, as this screws with plot dimensions)
    plt.subplots_adjust(left=options['subplots_adjust'][0], bottom=options['subplots_adjust'][1], right=options['subplots_adjust'][2], top=options['subplots_adjust'][3])


    # If limits for x- and y-axes is passed, sets these.
    if options['xlim']:
        ax.set_xlim(options['xlim'])

    if options['ylim']:
        ax.set_ylim(options['ylim'])


    # Add custom text
    if options['text']:
        plt.text(x=options['text'][1][0], y=options['text'][1][1], s=options['text'][0])



    if options['e_fermi']:
        if options['flip_xy']:
            ax.axhline(0, c='black', ls='dashed')
        else:
            ax.axvline(0, c='black', ls='dashed')

    if options['plot_kind'] == 'DOS' or options['plot_kind'] == 'PDOS':
        if options['dos_info']['spin_polarised']:
            if options['flip_xy']:
                ax.axvline(0, c='black')
            else:
                ax.axhline(0, c='black')
    elif options['plot_kind'] == 'COOP' or options['plot_kind'] == 'COHP':
        if options['flip_xy']:
            ax.axvline(0, c='black')
        else:
            ax.axhline(0, c='black')
    
    return fig, ax



def plot_coop(plot_data, options):
    ''' interactions = list with number of interaction (index + 1 of interactions list from read_coop)'''


    required_options = ['plot_kind', 'mode', 'up', 'down', 'collapse', 'interactions', 'flip_xy', 'fill', 'colours', 'palettes']

    default_options = {
        'plot_kind': 'COOP',
        'mode': 'individual',
        'fill': False,
        'up': True,
        'down': True,
        'collapse': False,
        'interactions': None,
        'palettes': [('qualitative', 'Dark2_8')],
        'colours': None,
        'flip_xy': False
    
    }

    options = update_options(options=options, required_options=required_options, default_options=default_options)


    fig, ax = prepare_plot(options=options)

    coopcar, coop_interactions = dft.io.read_coop(plot_data=plot_data, options=options)



    if not options['colours']:
        colour_cycle = generate_colours(palettes=options['palettes'])

        colours = []
        for interaction in range(len(coop_interactions)):
            colours.append(next(colour_cycle))

    else:
        colours = options['colours']


    # If interactions has been specified
    if options['interactions']:

        # Make interactions into a list of lists for correct looping below
        if type(options['interactions'][0]) != list:
            interactions_list = [options['interactions']]
        else:
            interactions_list = options['interactions']

        for ind, interactions in enumerate(interactions_list):

            # Determine which columns to plot if collapse is enabled
            if options['collapse']:
                to_plot = [2*(i-1)+3 for i in interactions]

                # Make sum column for plotting if sum mode is enabled
                if options['mode'] == 'sum':
                    coopcar["sum"] = coopcar.iloc[:, to_plot].sum(axis=1)
                    to_plot = ['sum']


                # Make mean column for plotting if mean mode is enabeld
                elif options['mode'] == 'avg' or options['mode'] == 'average' or options['mode'] == 'mean':
                    coopcar["mean"] = coopcar.iloc[:, to_plot].mean(axis=1)
                    to_plot = ['mean']                   
            
            # Determine which columns to plot if collapse is disabled and both up and down should be plotted
            elif options['up'] and options['down']:
                to_plot_up = [2*(i-1)+3 for i in interactions]
                to_plot_down = [2*(i-1)+5 +2*len(coop_interactions) for i in interactions]
                to_plot = to_plot_up + to_plot_down

                if options['mode'] == 'sum':
                    coopcar["sum_up"] = coopcar.iloc[:, to_plot_up].sum(axis=1)
                    coopcar["sum_down"] = coopcar.iloc[:, to_plot_down].sum(axis=1)
                    to_plot = ['sum_up', 'sum_down']

                elif options['mode'] == 'avg' or options['mode'] == 'average' or options['mode'] == 'mean':
                    coopcar["mean_up"] = coopcar.iloc[:, to_plot_up].mean(axis=1)
                    coopcar["mean_down"] = coopcar.iloc[:, to_plot_down].mean(axis=1)
                    to_plot = ['mean_up', 'mean_down']

            # Determine which columns to plot if collapse is disabled and only up should be plotted
            elif options['up']:
                to_plot = [2*(i-1)+3 for i in interactions]

                if options['mode'] == 'sum':
                    coopcar["sum_up"] = coopcar.iloc[:, to_plot].sum(axis=1)
                    to_plot = ['sum_up']

                elif options['mode'] == 'avg' or options['mode'] == 'average' or options['mode'] == 'mean':
                    coopcar["mean_up"] = coopcar.iloc[:, to_plot].mean(axis=1)
                    to_plot = ['mean_up']


             # Determine which columns to plot if collapse is disabled and only down should be plotted
            elif options['down']:
                to_plot = [2*(i-1)+5 +2*len(coop_interactions) for i in interactions]

                if options['mode'] == 'sum':
                    coopcar["sum_down"] = coopcar.iloc[:, to_plot].sum(axis=1)
                    to_plot = ['sum_down']

                elif options['mode'] == 'avg' or options['mode'] == 'average' or options['mode'] == 'mean':
                    coopcar["mean_down"] = coopcar.iloc[:, to_plot].mean(axis=1)
                    to_plot = ['mean_down']
                        

            # Plot all columns as decided above
            for j, column in enumerate(to_plot):
                if options['fill']:
                    ax.fill_between(coopcar["Energy"], coopcar[column], 0, where=coopcar[column]>0, color=colours[ind])
                    ax.fill_between(coopcar["Energy"], coopcar[column], 0, where=coopcar[column]<0, color=colours[ind+1])

                else:
                    if options['mode'] == "individual":
                        colour = colours[j]
                    else:
                        colour = colours[ind]

                    
                    if options['flip_xy']:
                        coopcar.plot(y='Energy', x=column, ax=ax, color=colour)
                    else:
                        coopcar.plot(x='Energy', y=column, ax=ax, color=colour)

    prettify_dos_plot(fig=fig, ax=ax, options=options)

    return coopcar



def prettify_coop_plot(fig, ax, energyunit='eV', dosunit='arb. u.', xlim=None, ylim=None, title=None, hide_ylabels=False, hide_xlabels=False, hide_yvals=False, hide_xvals=False, flip_xy=False, pad_bottom=None, scale=1, colours=None, atoms=None, pdos=False, width=None, height=None, e_fermi=False, adjust=False, legend=False, labels=None, label_colours=None, xpad=0, ypad=0):
    
    # Set sizes of ticks, labes etc.
    ticksize = 30*scale
    labelsize = 30*scale
    legendsize = 30*scale
    titlesize = 30*scale
    
    linewidth = 3*scale
    axeswidth = 3*scale
    majorticklength = 20*scale
    minorticklength = 10*scale
    
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    
    if flip_xy:
    
        # Set labels on x- and y-axes
        if not hide_ylabels:
            if ypad:
                ax.set_ylabel('Energy [{}]'.format(energyunit), size=labelsize, labelpad=ypad)
            else:
                ax.set_ylabel('Energy [{}]'.format(energyunit), size=labelsize)
        
        if pdos:
            if xpad:
                ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize, labelpad=xpad)
            else:
                ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize)

        else:
            if width >= 10:
                if xpad:
                    ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize, labelpad=xpad)
                else:
                    ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize)                    
            
            else:
                if xpad:
                    ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize, labelpad=xpad)
                else:
                    ax.set_xlabel('COOP [{}]'.format(dosunit), size=labelsize)                    

        ax.tick_params(axis='y', direction='in', which='major', right=True, length=majorticklength, width=linewidth)
        ax.tick_params(axis='y', direction='in', which='minor', right=True, length=minorticklength, width=linewidth)

        if hide_yvals:
            ax.tick_params(axis='y', labelleft=False)

        ax.tick_params(axis='x', direction='in', which='major', bottom=False, labelbottom=False)
    
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(.5))
                                      
    
    
    else:
        # Set labels on x- and y-axes
        if adjust:
            if xpad:
                ax.set_xlabel('E - E$_F$ [{}]'.format(energyunit), size=labelsize, labelpad=xpad)
            else:
                ax.set_xlabel('E - E$_F$ [{}]'.format(energyunit), size=labelsize)        


        else:
            if xpad:
                ax.set_xlabel('Energy [{}]'.format(energyunit), size=labelsize, labelpad=xpad)
            else:
                ax.set_xlabel('Energy [{}]'.format(energyunit), size=labelsize)                
        
        
        if height < 10:
            if ypad:
                ax.set_ylabel('COOP [{}]'.format(dosunit), size=labelsize, labelpad=ypad)
            else:
                ax.set_ylabel('COOP [{}]'.format(dosunit), size=labelsize)

        else:
            if ypad:
                ax.set_ylabel('Crystal orbital overlap population [{}]'.format(dosunit), size=labelsize, labelpad=ypad)
            else:
                ax.set_ylabel('Crystal orbital overlap population [{}]'.format(dosunit), size=labelsize)                

    
        ax.tick_params(axis='x', direction='in', which='major', bottom=True, top=True, length=majorticklength, width=linewidth)
        ax.tick_params(axis='x', direction='in', which='minor', bottom=True, top=True, length=minorticklength, width=linewidth)
    
 
        ax.tick_params(axis='y', which='major', direction='in', right=True, left=True, labelleft=True, length=majorticklength, width=linewidth)
        ax.tick_params(axis='y', which='minor', direction='in', right=True, left=True, length=minorticklength, width=linewidth)

        if hide_ylabels:
            ax.set_ylabel('')
        if hide_xlabels:
            ax.set_xlabel('')
        if hide_yvals:
            ax.tick_params(axis='y', which='both', labelleft=False)
        if hide_xvals:
            ax.tick_params(axis='x', which='both', labelbottom=False)
            

        if ylim:
            yspan = ylim[1] - ylim[0]
            yloc = np.round(yspan / 4, 2)
            
            ax.yaxis.set_major_locator(MultipleLocator(yloc))
            ax.yaxis.set_minor_locator(MultipleLocator(yloc/2))
    


        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.5)) 
    
    
    
    plt.xlim(xlim)
    plt.ylim(ylim)


    if title:
        ax.set_title(title, size=40)



    if legend:
        patches = []

        if label_colours:
            colours=label_colours

        for ind, label in enumerate(labels):
            patches.append(mpatches.Patch(color=colours[ind], label=label))
            
        fig.legend(handles=patches, loc='upper right', ncol=len(labels), bbox_to_anchor=(0.8, 0.45), fontsize=legendsize/1.25, frameon=False)

            #bbox_to_anchor=(1.20, 0.91)

        
        
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

    if e_fermi:
        if flip_xy:
            plt.axhline(0, lw=linewidth, c='black', ls='dashed')
        else:
            plt.axvline(0, lw=linewidth, c='black', ls='dashed')


    plt.axhline(0, lw=linewidth, c='black')
    
    return fig, ax




def get_unique_atoms(interactions):
    ''' Get all the unique atoms involved in the interactions from the COOP-calculation

    Input:
    interactions: list of interactions that comes as output from read_coop()

    Outut:
    unique_atoms: list of unique atoms in the interactions list'''

    unique_atoms = []

    for interaction in interactions:

        atoms = interaction.split('->')

        for atom in atoms:
            if atom not in unique_atoms:
                unique_atoms.append(atom)


    unique_atoms.sort()

    return unique_atoms


def get_interactions_involving(interactions, targets):
    ''' Get the indicies (+1) of all the interactions involving target. This list can be used as input to plot_coop(), as it is
    then formatted the way that function accepts these interactions.

    Input:
    interactions: list of interactions as output from read_coop()
    target: the particular atom that should be involved in the interactions contained in the output list

    Output:
    target_interactions: Indices (+1) of all the interactions involving target atom.'''

    target_interactions = []
    appended_interactions = []


    if type(targets) == list:
        for target in targets:
            for ind, interaction in enumerate(interactions):
                if target in interaction.split('->') and interaction not in appended_interactions:
                    target_interactions.append(ind+1)
                    appended_interactions.append(interaction)

    else:
        for ind, interaction in enumerate(interactions):
            if targets in interaction.split('->'):
                target_interactions.append(ind+1)


    return target_interactions




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



def generate_colours(options: dict):

    if not isinstance(options['palettes'], list):
        options['palettes'] = [options['palettes']]

    # Creates a list of all the colours that is passed in the colour_cycles argument. Then makes cyclic iterables of these. 
    colour_collection = []

    for palette in options['palettes']:
        mod = importlib.import_module("palettable.colorbrewer.%s" % palette[0])
        colour = getattr(mod, palette[1]).mpl_colors
        colour_collection = colour_collection + colour

    colour_cycle = itertools.cycle(colour_collection)


    return colour_cycle