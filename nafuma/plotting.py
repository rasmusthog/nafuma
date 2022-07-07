import nafuma.auxillary as aux

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import importlib
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import itertools


import numpy as np



def prepare_plot(options={}):
    ''' A general function to prepare a plot based on contents of options['rc_params'] and options['format_params'].
    
    rc_params is a dictionary with keyval-pairs corresponding to rcParams in matplotlib, to give the user full control over this. Please consult the matplotlib-documentation
    
    format_params will determine the size, aspect ratio, resolution etc. of the figure. Should be modified to conform with any requirements from a journal.'''

    if 'rc_params' in options.keys():
        rc_params = options['rc_params']
    else:
        rc_params = {}

    if 'format_params' in options.keys():
        format_params = options['format_params']
    else:
        format_params = {}
    
    required_format_params = ['single_column_width', 'double_column_width', 'column_type', 'width_ratio', 'aspect_ratio', 
    'width', 'height', 'compress_width', 'compress_height', 'upscaling_factor', 'dpi',
    'nrows', 'ncols', 'grid_ratio_height', 'grid_ratio_width']

    default_format_params = {
    'single_column_width': 8.3,
    'double_column_width': 17.1,
    'column_type': 'single',
    'width_ratio': '1:1',
    'aspect_ratio': '1:1',
    'width': None,
    'height': None,
    'compress_width': 1,
    'compress_height': 1,
    'upscaling_factor': 1.0,
    'dpi': 600,
    'nrows': 1,
    'ncols': 1,
    'grid_ratio_height': None,
    'grid_ratio_width': None
    }
    
    format_params = aux.update_options(format_params, required_format_params, default_format_params)


    # Reset run commands
    plt.rcdefaults()
    
    # Update run commands if any is passed (will pass an empty dictionary if not passed)
    update_rc_params(rc_params)
    
    if not format_params['width']:
        format_params['width'] = determine_width(format_params=format_params)
    
    if not format_params['height']:
        format_params['height'] = determine_height(format_params=format_params, width=format_params['width'])

    format_params['width'], format_params['height'] = scale_figure(format_params=format_params, width=format_params['width'], height=format_params['height'])
    
    if format_params['nrows'] == 1 and format_params['ncols'] == 1:
        fig, ax = plt.subplots(figsize=(format_params['width'], format_params['height']), dpi=format_params['dpi'])
        
        return fig, ax

    else:
        if not format_params['grid_ratio_height']:
            format_params['grid_ratio_height'] = [1 for i in range(format_params['nrows'])]

        if not format_params['grid_ratio_width']:
            format_params['grid-ratio_width'] = [1 for i in range(format_params['ncols'])]

        fig, axes = plt.subplots(nrows=format_params['nrows'], ncols=format_params['ncols'], figsize=(format_params['width'],format_params['height']), 
        gridspec_kw={'height_ratios': format_params['grid_ratio_height'], 'width_ratios': format_params['grid_ratio_width']}, 
        facecolor='w', dpi=format_params['dpi'])

        return fig, axes


def adjust_plot(fig, ax, options):
    ''' A general function to adjust plot according to contents of the options-dictionary '''
    
    required_options = [
    'plot_kind', 
    'xlabel', 'ylabel',
    'xunit', 'yunit',
    'hide_x_labels', 'hide_y_labels',
    'hide_x_ticklabels', 'hide_y_ticklabels',
    'hide_x_ticks', 'hide_y_ticks',
    'x_tick_locators', 'y_tick_locators',  
    'rotation_x_ticks', 'rotation_y_ticks',
    'xticks', 'yticks', 
    'xlim', 'ylim', 
    'title', 
    'legend', 'legend_position', 'legend_ncol',
    'subplots_adjust',
    'text']

    default_options = {
        'plot_kind': None, # defaults to None, but should be utilised when requiring special formatting for a particular plot 
        'xlabel': None, 'ylabel': None,
        'xunit': None, 'yunit': None,
        'hide_x_labels': False, 'hide_y_labels': False, # Whether the main labels on the x- and/or y-axes should be hidden
        'hide_x_ticklabels': False, 'hide_y_ticklabels': False, # Whether ticklabels on the x- and/or y-axes should be hidden
        'hide_x_ticks': False, 'hide_y_ticks': False, # Whether the ticks on the x- and/or y-axes should be hidden
        'x_tick_locators': None, 'y_tick_locators': None, # The major and minor tick locators for the x- and y-axes
        'rotation_x_ticks': 0, 'rotation_y_ticks': 0, # Degrees the x- and/or y-ticklabels should be rotated
        'xticks': None, 'yticks': None, # Custom definition of the xticks and yticks. This is not properly implemented now.  
        'xlim': None, 'ylim': None, # Limits to the x- and y-axes
        'title': None, # Title of the plot
        'legend': False, 'legend_position': ['lower center', (0.5, -0.1)], 'legend_ncol': 1, # Toggles on/off legend. Specifices legend position and the number of columns the legend should appear as.
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9], # Adjustment of the Axes-object within the Figure-object. Fraction of the Figure-object the left, bottom, right and top edges of the Axes-object will start.
        'text': None # Text to show in the plot. Should be a list where the first element is the string, and the second is a tuple with x- and y-coordinates. Could also be a list of lists to show more strings of text.
    }



    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    # Set labels on x- and y-axes
    if not options['hide_y_labels']:
        ax.set_ylabel(f'{options["ylabel"]} [{options["yunit"]}]')
    else:
        ax.set_ylabel('')
        
    if not options['hide_x_labels']:
        if not options['xunit']:
            ax.set_xlabel(f'{options["xlabel"]}')
        else:
            ax.set_xlabel(f'{options["xlabel"]} [{options["xunit"]}]')
    else:
        ax.set_xlabel('')


    # Set multiple locators
    if options['y_tick_locators']:
        ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
        ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))

    if options['x_tick_locators']:
        ax.xaxis.set_major_locator(MultipleLocator(options['x_tick_locators'][0]))
        ax.xaxis.set_minor_locator(MultipleLocator(options['x_tick_locators'][1]))

    
    # FIXME THIS NEEDS REWORK FOR IT TO FUNCTION PROPERLY!
    #if options['xticks']:
    #    ax.set_xticks(np.arange(plot_data['start'], plot_data['end']+1))
    #    ax.set_xticklabels(options['xticks'])
    # else:
    #     ax.set_xticks(np.arange(plot_data['start'], plot_data['end']+1))
    #     ax.set_xticklabels([x/2 for x in np.arange(plot_data['start'], plot_data['end']+1)])
        
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
        ax.set_title(options['title'], fontsize=plt.rcParams['font.size'])

     

    # Create legend

    if ax.get_legend():
        ax.get_legend().remove()

    
    if options['legend']:
        

        # Make palette and linestyles from original parameters
        if not options['colours']:
            colours = generate_colours(palettes=options['palettes'])
        else:
            colours = itertools.cycle(options['colours'])
        

        markers = itertools.cycle(options['markers'])
        
        # Create legend
        active_markers = []
        active_labels = []

        for label in options['labels']:


            # Discard next linestyle and colour if label is _
            if label == '_':
                _ = next(colours)
                _ = next(markers)

            else:
                active_markers.append(mlines.Line2D([], [], markeredgecolor=next(colours), color=(1, 1, 1, 0), marker=next(markers)))
                active_labels.append(label)

    

        ax.legend(active_markers, active_labels, frameon=False, loc=options['legend_position'][0], bbox_to_anchor=options['legend_position'][1], ncol=options['legend_ncol'])
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

        # If only a single element, put it into a list so the below for-loop works.
        if isinstance(options['text'][0], str):
            options['text'] = [options['text']]

        # Plot all passed texts
        for text in options['text']:
            plt.text(x=text[1][0], y=text[1][1], s=text[0])
    
    return fig, ax




def ipywidgets_update(func, data, options={}, **kwargs):
    ''' A general ipywidgets update function that can be passed to ipywidgets.interactive. To use this, you can run:

    import ipywidgets as widgets
    import beamtime.plotting as btp

    w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(my_func), plot_data=widgets.fixed(plot_data), options=widgets.fixed(options), key1=widget1, key2=widget2, key3=widget3)

    where key1, key2, key3 etc. are the values in the options-dictionary you want widget control of, and widget1, widget2, widget3 etc. are widgets to control these values, e.g. widgets.IntSlider(value=1, min=0, max=10)
    '''


    # Update the options-dictionary with the values from the widgets
    for key in kwargs:
        options[key] = kwargs[key]

    # Call the function with the plot_data and options-dictionaries
    func(data=data, options=options)



def determine_width(format_params):
    ''' '''
    
    conversion_cm_inch = 0.3937008 # cm to inch
    
    if format_params['column_type'] == 'single':
        column_width = format_params['single_column_width']
    elif format_params['column_type'] == 'double':
        column_width = format_params['double_column_width']
        
    column_width *= conversion_cm_inch
    
    
    width_ratio = [float(num) for num in format_params['width_ratio'].split(':')]

    
    width = column_width * width_ratio[0]/width_ratio[1]

    
    return width


def determine_height(format_params, width):
    
    aspect_ratio = [float(num) for num in format_params['aspect_ratio'].split(':')]
    
    height = width/(aspect_ratio[0] / aspect_ratio[1])
    
    return height


def scale_figure(format_params, width, height):
    width = width * format_params['upscaling_factor'] * format_params['compress_width']
    height = height * format_params['upscaling_factor'] * format_params['compress_height']

    return width, height



def update_rc_params(rc_params):
    ''' Update all passed run commands in matplotlib'''
    
    if rc_params:
        for key in rc_params.keys():
            plt.rcParams.update({key: rc_params[key]})


def generate_colours(palettes, kind=None):

    if kind == 'single':
        colour_cycle = itertools.cycle(palettes)

    else:
        # Creates a list of all the colours that is passed in the colour_cycles argument. Then makes cyclic iterables of these. 
        colour_collection = []
        for palette in palettes:
            mod = importlib.import_module("palettable.colorbrewer.%s" % palette[0])
            colour = getattr(mod, palette[1]).mpl_colors
            colour_collection = colour_collection + colour

        colour_cycle = itertools.cycle(colour_collection)


    return colour_cycle