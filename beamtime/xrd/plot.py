import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math

import beamtime.xrd as xrd

import beamtime.auxillary as aux
import beamtime.plotting as btp


def plot_diffractogram(plot_data, options={}):
    ''' Plots a diffractogram.
    
    Input:
    plot_data (dict): Must include path = string to diffractogram data, and plot_kind = (recx, beamline, image)'''

    # Update options
    required_options = ['x_vals', 'y_vals', 'ylabel', 'xlabel', 'xunit', 'yunit', 'scatter', 'plot_kind', 'rc_params', 'format_params']

    default_options = {
        'x_vals': '2th', 
        'y_vals': 'I',
        'ylabel': 'Intensity', 'xlabel': '2theta', 
        'xunit': 'deg', 'yunit': 'a.u.',
        'scatter': False,
        'plot_kind': None,
        'rc_params': {},
        'format_params': {}
        }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    # Prepare plot, and read and process data
    fig, ax = btp.prepare_plot(options=options)
    diffractogram = xrd.io.read_data(path=plot_data['path'], kind=plot_data['plot_kind'], options=options)

    if options['scatter']:
        ax.scatter(x= diffractogram[options['x_vals']], y = diffractogram[options['y_vals']])
        
        #diffractogram.plot(x=options['x_vals'], y=options['y_vals'], ax=ax, kind='scatter')
    
    else:
        diffractogram.plot(x=options['x_vals'], y=options['y_vals'], ax=ax)



    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)


    return diffractogram, fig, ax



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