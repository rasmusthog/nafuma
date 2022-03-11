import beamtime.auxillary as aux
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
import importlib
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
from cycler import cycler
import itertools



def prepare_plot(options={}):
    ''' Prepares plot based on contents of options['rc_params'] and options['format_params'].
    
    rc_params is a dictionary with keyval-pairs corresponding to rcParams in matplotlib
    
    format_params will determine the size and aspect ratios of '''

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
    
    options = aux.update_options(format_params, required_options, default_options)


    # Reset run commands
    plt.rcdefaults()
    
    # Update run commands if any is passed (will pass an empty dictionary if not passed)
    update_rc_params(rc_params)
    
    width = determine_width(options)
    height = determine_height(options, width)
    width, height = scale_figure(options=options, width=width, height=height)
    
    fig, ax = plt.subplots(figsize=(width, height), dpi=options['dpi'])
    
    return fig, ax


def prettify__plot(fig, ax, plot_data, options):
    
    required_options = ['plot_kind', 'hide_x_labels', 'hide_y_labels',  'rotation_x_ticks', 'rotation_y_ticks', 'xlabel', 'ylabel', 'yunit', 'xlim', 'ylim', 'x_tick_locators', 'y_tick_locators', 'xticks', 'hide_x_ticks', 'hide_y_ticks', 'hide_x_ticklabels', 'hide_y_ticklabels',
                        'colours', 'palettes',  'title', 'legend', 'legend_position', 'subplots_adjust', 'text', 'legend_ncol']

    default_options = {
        'plot_kind': None, # defaults to None, but should be utilised when 
        'hide_x_labels': False, # Whether x labels should be hidden
        'hide_x_ticklabels': False,
        'hide_x_ticks': False,
        'rotation_x_ticks': 0,
        'hide_y_labels': False, # whether y labels should be hidden
        'hide_y_ticklabels': False,
        'hide_y_ticks': False,
        'rotation_y_ticks': 0,
        'xlabel': r'$x$ in Na$_{5-x}$FeO$_{4-\delta}$',
        'ylabel': 'Formation energy',
        'yunit': r'eV', # The unit of the y-values in the curve and bar plots
        'xlim': None,
        'ylim': None,
        'x_tick_locators': [.5, .25], # Major and minor tick locators
        'y_tick_locators': [.5, .25],
        'xticks': None,
        'labels': None,
        'colours': None,
        'palettes': [('qualitative', 'Dark2_8'), ('qualitative', 'Paired_12')],
        'title': None,
        'legend': False,
        'legend_position': ['lower center', (0.5, -0.1)], # the position of the legend passed as arguments to loc and bbox_to_anchor respectively 
        'legend_ncol': 1,
        'subplots_adjust': [0.1, 0.1, 0.9, 0.9],
        'text': None
    }



    options = update_options(options=options, required_options=required_options, default_options=default_options)

    # Set labels on x- and y-axes
    if not options['hide_y_labels']:
        ax.set_ylabel(f'{options["ylabel"]} [{options["yunit"]}]')
    else:
        ax.set_ylabel('')
        
    if not options['hide_x_labels']:
        ax.set_xlabel(f'{options["xlabel"]}')
    else:
        ax.set_xlabel('')


    # Set multiple locators
    ax.yaxis.set_major_locator(MultipleLocator(options['y_tick_locators'][0]))
    ax.yaxis.set_minor_locator(MultipleLocator(options['y_tick_locators'][1]))

    ax.xaxis.set_major_locator(MultipleLocator(options['x_tick_locators'][0]))
    ax.xaxis.set_minor_locator(MultipleLocator(options['x_tick_locators'][1]))

    if options['xticks']:
        ax.set_xticks(np.arange(plot_data['start'], plot_data['end']+1))
        ax.set_xticklabels(options['xticks'])
    else:
        ax.set_xticks(np.arange(plot_data['start'], plot_data['end']+1))
        ax.set_xticklabels([x/2 for x in np.arange(plot_data['start'], plot_data['end']+1)])
        
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
        plt.text(x=options['text'][1][0], y=options['text'][1][1], s=options['text'][0])
    
    return fig, ax




def ipywidgets_update(func, plot_data, options={}, **kwargs):

    for key in kwargs:
        options[key] = kwargs[key]

    func(plot_data=plot_data, options=options)




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



def update_rc_params(rc_params):
    ''' Update all passed run commands in matplotlib'''
    
    if rc_params:
        for key in rc_params.keys():
            plt.rcParams.update({key: rc_params[key]})


def generate_colours(palettes):

    # Creates a list of all the colours that is passed in the colour_cycles argument. Then makes cyclic iterables of these. 
    colour_collection = []
    for palette in palettes:
        mod = importlib.import_module("palettable.colorbrewer.%s" % palette[0])
        colour = getattr(mod, palette[1]).mpl_colors
        colour_collection = colour_collection + colour

    colour_cycle = itertools.cycle(colour_collection)


    return colour_cycle