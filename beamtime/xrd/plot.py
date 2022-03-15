import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math

import ipywidgets as widgets

import beamtime.xrd as xrd
import beamtime.auxillary as aux
import beamtime.plotting as btp


def plot_diffractogram(plot_data, options={}):
    ''' Plots a diffractogram.
    
    Input:
    plot_data (dict): Must include path = string to diffractogram data, and plot_kind = (recx, beamline, image)'''

    # Update options
    required_options = ['x_vals', 'y_vals', 'ylabel', 'xlabel', 'xunit', 'yunit', 'line', 'scatter', 
    'reflections_plot', 'reflections_indices', 'reflections_data', 'plot_kind', 'palettes', 'interactive', 'rc_params', 'format_params']

    default_options = {
        'x_vals': '2th', 
        'y_vals': 'I',
        'ylabel': 'Intensity', 'xlabel': '2theta', 
        'xunit': 'deg', 'yunit': 'a.u.',
        'xlim': None, 'ylim': None, 
        'line': True, # whether or not to plot diffractogram as a line plot
        'scatter': False, # whether or not to plot individual data points
        'reflections_plot': False, # whether to plot reflections as a plot
        'reflections_indices': False, # whether to plot the reflection indices
        'reflections_data': None, # Should be passed as a list of dictionaries on the form {path: rel_path, reflection_indices: number of indices, colour: [r,g,b], min_alpha: 0-1]
        'plot_kind': None,
        'palettes': [('qualitative', 'Dark2_8')],
        'interactive': False,
        'interactive_session_active': False,
        'rc_params': {},
        'format_params': {},
        }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if options['interactive']:
        options['interactive'] = False
        options['interactive_session_active'] = True
        plot_diffractogram_interactive(plot_data=plot_data, options=options)
        return
    
    
    if options['reflections_data']:
        if not isinstance(options['reflections_data'], list):
            options['reflections_data'] = [options['reflections_data']]

    # Determine number of subplots and height ratios between them
    if len(options['reflections_data']) >= 1:
        options = determine_grid_layout(options=options)

        print(options['format_params']['nrows'])


    # Prepare plot, and read and process data

    fig, ax = btp.prepare_plot(options=options)
    

    # Assign the correct axes
    if options['reflections_plot'] or options['reflections_indices']:
        
        if options['reflections_indices']:
            indices_ax = ax[0]

            if options['reflections_plot']:
                ref_axes = [axx for axx in ax[range(1,len(options['reflections_data'])+1)]]

        else:
            ref_axes = [axx for axx in ax[range(0,len(options['reflections_data']))]]

        ax = ax[-1]

    colours = btp.generate_colours(options['palettes'])


    diffractogram = xrd.io.read_data(path=plot_data['path'], kind=plot_data['plot_kind'], options=options)


    if options['line']:
        diffractogram.plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=next(colours), zorder=1)
    
    if options['scatter']:
        ax.scatter(x=diffractogram[options['x_vals']], y = diffractogram[options['y_vals']], c=[(1,1,1,0)], edgecolors=[next(colours)], linewidths=plt.rcParams['lines.markeredgewidth'], zorder=2) #, edgecolors=np.array([next(colours)]))



    if not options['xlim']:
        options['xlim'] = [diffractogram[options['x_vals']].min(), diffractogram[options['x_vals']].max()]

    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)

    

    # Make the reflection plots
    if options['reflections_plot'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        
        for reference, axis in zip(options['reflections_data'], ref_axes):
            plot_reflection_table(plot_data=reference, ax=axis, options=options)

    # Print the reflection indices
    if options['reflections_indices'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        for reference in options['reflections_data']:
            plot_reflection_indices(plot_data=reference, ax=indices_ax, options=options)


    return diffractogram, fig, ax


def determine_grid_layout(options):


    nrows = 1 if not options['reflections_indices'] else 2

    if options['reflections_plot']:
        for reference in options['reflections_data']:
            nrows += 1

    options['format_params']['nrows'] = nrows
    options['format_params']['grid_ratio_height'] = [1 for i in range(nrows-1)]+[10]

    return options



def plot_diffractogram_interactive(plot_data, options):

    if options['reflections_data']:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), plot_data=widgets.fixed(plot_data), options=widgets.fixed(options), 
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True), 
        reflections_plot=widgets.ToggleButton(value=True),
        reflections_indices=widgets.ToggleButton(value=False))
    
    else:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), plot_data=widgets.fixed(plot_data), options=widgets.fixed(options), 
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True))
    
    
    display(w)


def plot_reflection_indices(plot_data, ax, options={}):
    ''' Print reflection indices from output generated by VESTA.
    
    Required contents of plot_data:
    path (str): relative path to reflection table file'''

    required_options = ['reflection_indices', 'text_colour', 'hide_indices']

    default_options = {
        'reflection_indices': 3, # Number of reflection indices to plot, from highest intensity and working its way down
        'text_colour': 'black',
        'hide_indices': False
    }

    plot_data = update_options(options=plot_data, required_options=required_options, default_options=default_options)

    if not plot_data['hide_indices']:
        reflection_table = xrd.io.load_reflection_table(plot_data['path'])
        
        if plot_data['reflection_indices'] > 0:
            reflection_indices = reflection_table.nlargest(options['reflection_indices'], 'I')

        
        for i in range(plot_data['reflection_indices']):
            ax.text(s=f'({reflection_indices["h"].iloc[i]} {reflection_indices["k"].iloc[i]} {reflection_indices["l"].iloc[i]})', x=reflection_indices['2th'].iloc[i], y=0, fontsize=2.5, rotation=90, va='bottom', ha='center', c=plot_data['text_colour'])    

        
        if options['xlim']:
            ax.set_xlim(options['xlim'])
        
        ax.axis('off')


    return
    
def plot_reflection_table(plot_data, ax=None, options={}):
    ''' Plots a reflection table from output generated by VESTA.
    
    Required contents of plot_data:
    path (str): relative path to reflection table file'''

    required_options = ['reflection_indices', 'reflections_colour', 'min_alpha', 'format_params', 'rc_params', 'label']

    default_options = {
        'reflection_indices': 0, # Number of indices to print
        'reflections_colour': [0,0,0],
        'min_alpha': 0,
        'format_params': {},
        'rc_params': {},
        'label': None
    }

    if 'colour' in plot_data.keys():
        options['reflections_colour'] = plot_data['colour']
    if 'min_alpha' in plot_data.keys():
        options['min_alpha'] = plot_data['min_alpha']
    if 'reflection_indices' in plot_data.keys():
        options['reflection_indices'] = plot_data['reflection_indices']
    if 'label' in plot_data.keys():
        options['label'] = plot_data['label']


    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    
    if not ax:
        _, ax = btp.prepare_plot(options)

    reflection_table = xrd.io.load_reflection_table(plot_data['path'])

    reflections, intensities  = reflection_table['2th'], reflection_table['I']


    
    colours = []

    for ref, intensity in zip(reflections, intensities):

        colour = list(options['reflections_colour'])
        rel_intensity = (intensity / intensities.max())*(1-options['min_alpha']) + options['min_alpha']
        colour.append(rel_intensity)
        colours.append(colour)

    

    
    ax.vlines(x=reflections, ymin=-1, ymax=1, colors=colours, lw=0.5)
    ax.set_ylim([-0.5,0.5])


    ax.tick_params(which='both', bottom=False, labelbottom=False, right=False, labelright=False, left=False, labelleft=False, top=False, labeltop=False)
    
    if options['xlim']:
        ax.set_xlim(options['xlim'])


    if options['label']:
        xlim_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim_avg = (ax.get_ylim()[0]+ax.get_ylim()[1])/2

        ax.text(s=plot_data['label'], x=(ax.get_xlim()[0]-0.01*xlim_range), y=ylim_avg, ha = 'right', va = 'center')




    


def prepare_diffractogram_plot(options=None):
    # First take care of the options for plotting - set any values not specified to the default values
    required_options = ['columns', 'width', 'height', 'format', 'dpi',  'facecolor']
    default_options = {'columns': 1, 'width': 14, 'format': 'golden_ratio', 'dpi': None, 'facecolor': 'w'}


    # Define the required sizes
    required_sizes = ['lines', 'axes']




    # If none are set at all, just pass the default_options
    if not options:
        options = default_options
        options['height'] = options['width'] * (math.sqrt(5) - 1) / 2
        options['figsize'] = (options['width'], options['height'])

    
    # Define default sizes
    default_sizes = {
            'lines': 3*options['columns'],
            'axes': 3*options['columns']
    }

    # Initialise dictionary if it doesn't exist
    if not 'sizes' in options.keys():
        options['sizes'] = {}


    # Update dictionary with default values where none is supplied
    for size in required_sizes:
        if size not in options['sizes']:
            options['sizes'][size] = default_sizes[size]


    # If options is passed, go through to fill out the rest. 
    else:
        # Start by setting the width:
        if 'width' not in options.keys():
            options['width'] = default_options['width']

        # Then set height - check options for format. If not given, set the height to the width scaled by the golden ratio - if the format is square, set the same. This should possibly allow for the tweaking of custom ratios later.
        if 'height' not in options.keys():
            if 'format' not in options.keys():
                options['height'] = options['width'] * (math.sqrt(5) - 1) / 2
            elif options['format'] == 'square':
                options['height'] = options['width']

        options['figsize'] = (options['width'], options['height'])

        # After height and width are set, go through the rest of the options to make sure that all the required options are filled
        for option in required_options:
            if option not in options.keys():
                options[option] = default_options[option]

    fig, ax = plt.subplots(figsize=(options['figsize']), dpi=options['dpi'], facecolor=options['facecolor'])


    plt.rc('lines', linewidth=options['sizes']['lines'])
    plt.rc('axes', linewidth=options['sizes']['axes'])

    return fig, ax

def prettify_diffractogram_plot(fig, ax, options=None):

    ##################################################################
    ######################### UPDATE OPTIONS #########################
    ##################################################################

    # Define the required options
    required_options = [
    'columns', 
    'xticks', 'yticks', 
    'units',
    'show_major_ticks', 'show_minor_ticks', 'show_tick_labels',
    'xlim', 'ylim',
    'hide_x_axis', 'hide_y_axis', 
    'positions', 
    'xlabel', 'ylabel', 
    'sizes',
    'title'
    ]


    # Define the default options
    default_options = {
        'columns': 1,
        'xticks': [10, 5], 'yticks': [10000, 5000],
        'units': {'2th': '$^o$', 'I': 'arb. u.'},
        'show_major_ticks': [True, False, True, False], 'show_minor_ticks': [True, False, True, False], 'show_tick_labels': [True, False, False, False],
        'xlim': None,'ylim': None,
        'hide_x_axis': False, 'hide_y_axis': False,
        'positions': {'xaxis': 'bottom', 'yaxis': 'left'},
        'xlabel': None, 'ylabel': None,
        'sizes': None,
        'title': None 	
    }

    options = update_options(options, required_options, default_options)



    ##################################################################
    ########################## DEFINE SIZES ##########################
    ##################################################################

    # Define the required sizes
    required_sizes = [
        'labels', 
        'legend', 
        'title', 
        'line', 'axes', 
        'tick_labels', 
        'major_ticks', 'minor_ticks']



    # Define default sizes
    default_sizes = {
            'labels': 30*options['columns'],
            'legend': 30*options['columns'],
            'title': 30*options['columns'],
            'line': 3*options['columns'],
            'axes': 3*options['columns'],
            'tick_labels': 30*options['columns'],
            'major_ticks': 20*options['columns'],
            'minor_ticks': 10*options['columns']	
    }

    # Initialise dictionary if it doesn't exist
    if not options['sizes']:
        options['sizes'] = {}


    # Update dictionary with default values where none is supplied
    for size in required_sizes:
        if size not in options['sizes']:
            options['sizes'][size] = default_sizes[size]


    ##################################################################
    ########################## AXIS LABELS ###########################
    ##################################################################


    if not options['xlabel']:
        options['xlabel'] = prettify_labels(options['x_vals']) + ' [{}]'.format(options['units'][options['x_vals']])

    else:
        options['xlabel'] = options['xlabel'] + ' [{}]'.format(options['units'][options['x_vals']])


    if not options['ylabel']:
        options['ylabel'] = prettify_labels(options['y_vals']) + ' [{}]'.format(options['units'][options['y_vals']])

    else:
        options['ylabel'] = options['ylabel'] + ' [{}]'.format(options['units'][options['y_vals']])

    ax.set_xlabel(options['xlabel'], size=options['sizes']['labels'])
    ax.set_ylabel(options['ylabel'], size=options['sizes']['labels'])

    ##################################################################
    ###################### TICK MARKS & LABELS #######################
    ##################################################################

    ax.tick_params(
        direction='in', 
        which='major', 
        bottom=options['show_major_ticks'][0], labelbottom=options['show_tick_labels'][0],
        left=options['show_major_ticks'][1], labelleft=options['show_tick_labels'][1],
        top=options['show_major_ticks'][2], labeltop=options['show_tick_labels'][2],
        right=options['show_major_ticks'][3], labelright=options['show_tick_labels'][3],
        length=options['sizes']['major_ticks'],
        width=options['sizes']['axes'])


    ax.tick_params(direction='in', which='minor', bottom=options['show_minor_ticks'][0], left=options['show_minor_ticks'][1], top=options['show_minor_ticks'][2], right=options['show_minor_ticks'][3], length=options['sizes']['minor_ticks'], width=options['sizes']['axes'])



    if options['positions']['yaxis'] == 'right':
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()


    if options['hide_x_axis']:
        ax.axes.xaxis.set_visible(False)

    if options['hide_y_axis']:
        ax.axes.yaxis.set_visible(False)



    # Otherwise apply user input
    if options['xticks']:
        major_xtick = options['xticks'][0]
        minor_xtick = options['xticks'][1]


    if options['yticks']:

        major_ytick = options['yticks'][0]
        minor_ytick = options['yticks'][1]


    # Apply values
    ax.xaxis.set_major_locator(MultipleLocator(major_xtick))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_xtick))

    ax.yaxis.set_major_locator(MultipleLocator(major_ytick))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_ytick))




    # SET FONTSIZE OF TICK LABELS

    plt.xticks(fontsize=options['sizes']['tick_labels'])
    plt.yticks(fontsize=options['sizes']['tick_labels'])

    ##################################################################
    ########################## AXES LIMITS ###########################
    ##################################################################

    if options['xlim']:
        plt.xlim(options['xlim'])

    if options['ylim']:
        plt.ylim(options['ylim'])

    ##################################################################
    ############################# TITLE ##############################
    ##################################################################

    if options['title']:
        ax.set_title(options['title'], size=options['sizes']['title'])

    ##################################################################
    ############################# LEGEND #############################
    ##################################################################

    if ax.get_legend():
        ax.get_legend().remove()

    return fig, ax



def prettify_labels(label):
	
	labels_dict = {
		'2th': '2$\\theta$',
        'I': 'Intensity'	
        }

	return labels_dict[label]

def plot_diffractograms(paths, kind, options=None):


    fig, ax = prepare_diffractogram_plot(options=options)

    diffractograms = []

    for path in paths:
        diffractogram = xrd.io.read_data(path=path, kind=kind, options=options)
        diffractograms.append(diffractogram)


    required_options = ['type', 'xvals', 'yvals', 'x_offset', 'y_offset', 'normalise', 'normalise_around', 'reverse_order']
    default_options = {
        'type': 'stacked',
        'xvals': '2th',
        'yvals': 'I',
        'x_offset': 0,
        'y_offset': 0.2,
        'normalise': True,
        'normalise_around': None,
        'reverse_order': False
    }


    # If reverse_order is enabled, reverse the order
    if options['reverse_order']:
        diffractograms = reverse_diffractograms(diffractograms)


    # If normalise is enbaled, normalise all the diffractograms
    if options['normalise']:
        if not options['normalise_around']:
            for diffractogram in diffractograms:
                diffractogram["I"] = diffractogram["I"]/diffractogram["I"].max()
            else:
                diffractogram["I"] = diffractogram["I"]/diffractogram["I"].loc[(diffractogram['2th'] > options['normalise_around'][0]) & (diffractogram['2th'] < options['normalise_around'][1])].max()

    
    if options['type'] == 'stacked':
        for diffractogram in diffractograms:
            diffractogram.plot(x=options['xvals'], y=options['yvals'], ax=ax)


    fig, ax = prettify_diffractogram_plot(fig=fig, ax=ax, options=options)


    return diffractogram, fig, ax


def reverse_diffractograms(diffractograms):

    rev_diffractograms = []

    for i in len(diffractograms):
        rev_diffractograms.append(diffractograms.pop())

    return rev_diffractograms

#def plot_heatmap():





def update_options(options, required_options, default_options):

	if not options:
		options = default_options
	
	else:
		for option in required_options:
			if option not in options.keys():
				options[option] = default_options[option]

	return options