import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import pandas as pd
import numpy as np
import math
import os
import shutil

import ipywidgets as widgets
from IPython.display import display
from PIL import Image

import nafuma.xrd as xrd
import nafuma.auxillary as aux
import nafuma.plotting as btp

def plot_diffractogram(data, options={}):
    ''' Plots a diffractogram.
    
    Input:
    data (dict): Must include path = string to diffractogram data, and plot_kind = (recx, beamline, image)'''

    # Update options
    required_options = ['x_vals', 'y_vals', 'ylabel', 'xlabel', 'xunit', 'yunit', 'line', 'scatter', 'xlim', 'ylim', 'normalise', 'offset', 'offset_x', 'offset_y', 'offset_change',
    'reflections_plot', 'reflections_indices', 'reflections_data', 'heatmap', 'cmap', 'plot_kind', 'palettes', 'highlight', 'highlight_colours', 'interactive', 'rc_params', 'format_params', 'interactive_session_active', 'plot_diff']

    default_options = {
        'x_vals': '2th', 'y_vals': 'I',
        'xlabel': '2$\\theta$', 'ylabel': None, 
        'xunit': '$^{\circ}$', 'yunit': None,
        'xlim': None, 'ylim': None, 
        'normalise': True,
        'offset': True,
        'offset_x': 0,
        'offset_y': 1,
        'offset_change': False,
        'line': True, # whether or not to plot diffractogram as a line plot
        'scatter': False, # whether or not to plot individual data points
        'reflections_plot': False, # whether to plot reflections as a plot
        'reflections_indices': False, # whether to plot the reflection indices
        'reflections_data': None, # Should be passed as a list of dictionaries on the form {path: rel_path, reflection_indices: number of indices, colour: [r,g,b], min_alpha: 0-1]
        'heatmap': False,
        'cmap': 'viridis',
        'plot_kind': None,
        'palettes': [('qualitative', 'Dark2_8')],
        'highlight': None,
        'highlight_colours': ['red'],
        'interactive': False,
        'interactive_session_active': False,
        'rc_params': {},
        'format_params': {},
        'plot_diff': False,
        }

    if 'offset_y' not in options.keys():
        if len(data['path']) > 10:
            default_options['offset_y'] = 0.05

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    #options['current_offset_y'] = options['offset_y']

    # Convert data['path'] to list to allow iteration over this to accommodate both single and multiple diffractograms
    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

    
    ############################################################################################################################################################
    ##### LOADING DATA #########################################################################################################################################
    ############################################################################################################################################################
    
    # Check if there is some data stored already, load in data if not. This speeds up replotting in interactive mode.
    if not 'diffractogram' in data.keys():

        # This is to set the default values of the diffractogram y-label and -unit so that the actual yunit and ylable can switch back and forth between these and the heatmap values
        options['diff.yunit'] = 'a.u.'
        options['diff.ylabel'] = 'Intensity'

        # Initialise empty list for diffractograms and wavelengths. If wavelength is not manually passed it should be automatically gathered from the .xy-file
        data['diffractogram'] = [None for _ in range(len(data['path']))]

        if 'wavelength' not in data.keys():
            data['wavelength'] = [None for _ in range(len(data['path']))]
        else:
            # If only a single value is passed it should be set to be the same for all diffractograms passed
            if not isinstance(data['wavelength'], list):
                data['wavelength'] = [data['wavelength'] for _ in range(len(data['path']))] 

        
        
        
        # LOAD DIFFRACTOGRAMS
        
        if 'htxrd' in data.keys() and data['htxrd']:
            data['diffractogram'], data['wavelength'] = xrd.io.read_htxrd(data=data, options=options, index=0)
        
        else:
            for index in range(len(data['path'])):
                diffractogram, wavelength = xrd.io.read_data(data=data, options=options, index=index)


                data['diffractogram'][index] = diffractogram
                data['wavelength'][index] = wavelength


                # FIXME This is a quick fix as the image is not reloaded when passing multiple beamline datasets. Should probably be handled in io?
                data['image'] = None

        
        # Sets the xlim if this has not been specified
        if not options['xlim']:
            options['xlim'] = [data['diffractogram'][0][options['x_vals']].min(), data['diffractogram'][0][options['x_vals']].max()]
            
        # GENERATE HEATMAP DATA
        data['heatmap'], data['heatmap_xticks'], data['heatmap_xticklabels'], data['heatmap_yticks'], data['heatmap_yticklabels'] = generate_heatmap(data=data, options=options)
        options['heatmap_loaded'] = True

        if options['heatmap']:
            xlim_start_frac, xlim_end_frac = options['xlim'][0] / data['diffractogram'][0][options['x_vals']].max(), options['xlim'][1] / data['diffractogram'][0][options['x_vals']].max()
            options['xlim'] = [options['heatmap_xlim'][0]*xlim_start_frac, options['heatmap_xlim'][1]*xlim_end_frac]

        if options['heatmap_reverse']:
            data['heatmap'] = data['heatmap'].iloc[::-1]
            data['heatmap_yticklabels'] = data['heatmap_yticklabels'][::-1]

    # If data was already loaded, only do a check to see if the data is in a list or not, and if not, put it in one. This is because it will be looped over later.
    else:
        if not isinstance(data['diffractogram'], list):
            data['diffractogram'] = [data['diffractogram']]
            data['wavelength'] = [data['wavelength']]

    ############################################################################################################################################################
    ##### INTERACTIVE SESSION ##################################################################################################################################
    ############################################################################################################################################################



    # START INTERACTIVE SESSION
    # Start inteactive session with ipywidgets. Disables options['interactive'] in order for the interactive loop to not recursively start new interactive sessions
    if options['interactive']:
        options['interactive'] = False
        options['interactive_session_active'] = True
        plot_diffractogram_interactive(data=data, options=options)
        return
    
    
    # If interactive mode is already enabled, update the offsets. 
    if options['interactive_session_active']:
        if options['offset']:
            if (options['offset_x'] != options['current_offset_x']) or (options['offset_y'] != options['current_offset_y']):
                for i, (diff, wl) in enumerate(zip(data['diffractogram'], data['wavelength'])):
                    xrd.io.apply_offset(diff, wl, i, options)

    


    ############################################################################################################################################################
    ##### PREPARE THE PLOT AND COLOURS #########################################################################################################################
    ############################################################################################################################################################
    
    # CREATE AND ASSIGN AXES

    # Makes a list out of reflections_data if it only passed as a dict, as it will be looped through later
    if options['reflections_data']:
        if not isinstance(options['reflections_data'], list):
            options['reflections_data'] = [options['reflections_data']]
    
    
    # Determine the grid layout based on how many sets of reflections data has been passed
    if options['reflections_data'] and len(options['reflections_data']) >= 1:
        options = determine_grid_layout(options=options)

    # Create the Figure and Axes objects
    fig, ax = btp.prepare_plot(options=options)

    # Assign the correct axes to the indicies, reflections and figure itself
    if options['reflections_plot'] or options['reflections_indices']:
        
        if options['reflections_indices']:
            indices_ax = ax[0]

            if options['reflections_plot']:
                ref_axes = [axx for axx in ax[range(1,len(options['reflections_data'])+1)]]

        else:
            ref_axes = [axx for axx in ax[range(0,len(options['reflections_data']))]]

        ax = ax[-1]

    
    # GENERATE COLOURS
    
    # Limit for when it is assumed that each diffractogram should have its own colour - after 8, the default colour palette is used up and starts a new.
    # FIXME Should probably allow for more than 8 if wanted - not a priority now
    if len(data['path']) <= 8:
        if 'colours' in options.keys():
            colours = btp.generate_colours(options['colours'], kind='single')

        else:
            colours = btp.generate_colours(options['palettes'])
    

    # Generates the colours of a list of scans to highlight is passed. options['highlight'] and options['highlight_colour'] must be of equal length. Entries in highlight can either be a list or a single number,
    # if the latter it will be turned into a list with the same number as element 1 and 2. 
    elif options['highlight']:
        # Make sure that options['highlight'] is a list
        if not isinstance(options['highlight'], list):
            options['highlight'] = [[options['highlight'], options['highlight']]]
        
        # Make sure that options['highlight_colours] is a list
        if not isinstance(options['highlight_colours'], list):
            options['highlight_colours'] = [options['highlight_colours']]

        colours = []
    
        # Loop through each scan - assign the correct colour to each of the scan intervals in options['highlight']
        for i in range(len(data['path'])):
            assigned = False
            for j, highlight in enumerate(options['highlight']):
                
                # If one of the elements in options['highlight'] is a single number (i.e. only one scan should be highlighted), this is converted into the suitable format to be handled below
                if not isinstance(highlight, list):
                    highlight = [highlight, highlight]

                # Assigns the j-th colour if scan number (i) is within the j-th highlight-interval
                if i >= highlight[0] and i <= highlight[1]:
                    colours.append(options['highlight_colours'][j])
                    assigned = True
            
            # Only assign black to i if not already been given a colour
            if not assigned:
                colours.append('black')

            # Reset the 'assigned' value for the next iteration
            assigned = False

        # Make a itertools cycle out of the colours
        colours = btp.generate_colours(colours, kind='single')


    # If there are many scans and no highlight-options have been passed, all scans will be black
    else:
        colours = btp.generate_colours(['black'], kind='single')



    ############################################################################################################################################################
    ##### PLOT THE DATA ########################################################################################################################################
    ############################################################################################################################################################


    # PLOT HEATMAP
    if options['heatmap']:

        # Add locators for y-axis - otherwise it will tend to break (too many ticks) when switching between diffractograms and heatmap in interactive mode. These values will be updated later anyway, and is only 
        # to allow the initial call to Seaborn to have values that are sensible.
        # FIXME A more elegant solution to this?
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(50))

        # Call Seaborn to plot the data
        sns.heatmap(data['heatmap'], cmap=options['cmap'], cbar=False, ax=ax)
     
        
        # Set the ticks and ticklabels to match the data point number with 2th values 
        ax.set_xticks(data['heatmap_xticks'][options['x_vals']])
        ax.set_xticklabels(data['heatmap_xticklabels'][options['x_vals']])
        ax.set_yticks(data['heatmap_yticks'])
        ax.set_yticklabels(data['heatmap_yticklabels'])

        # Set the labels to the relevant values for heatmap plot
        if not options['ylabel'] or options['ylabel'] == options['diff.ylabel']:
            options['ylabel'] = options['heatmap.ylabel']
        if not options['yunit'] or options['yunit'] == options['diff.yunit']:
            options['yunit'] = options['heatmap.yunit']
        
        

        ax.tick_params(axis='x', which='minor', bottom=False, top=False)
        ax.tick_params(axis='y', which='minor', left=False, right=False)

        options['hide_y_ticklabels'] = False
        options['hide_y_ticks'] = False


        # Toggle on the frame around the heatmap - this makes it look better together with axes ticks
        for _, spine in ax.spines.items():
            spine.set_visible(True)


        if options['highlight']:
            for i, highlight in enumerate(options['highlight']):
                if i < len(options['highlight']) or len(options['highlight']) == 1: 
                    ax.axhline(y=highlight[1], c=options['highlight_colours'][i], ls='--', lw=0.5)


    # PLOT DIFFRACTOGRAM
    else:
        for diffractogram in data['diffractogram']:

            # Plot data as line plot
            if options['line']:
                diffractogram.plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=next(colours), zorder=1)
        
            # Plot data as scatter plot
            if options['scatter']:
                ax.scatter(x=diffractogram[options['x_vals']], y = diffractogram[options['y_vals']], c=[(1,1,1,0)], edgecolors=[next(colours)], linewidths=plt.rcParams['lines.markeredgewidth'], zorder=2) #, edgecolors=np.array([next(colours)]))


        # Set the labels to the relevant values for diffractogram plot
        if not options['ylabel'] or options['ylabel'] == options['heatmap.ylabel']:
            options['ylabel'] = options['diff.ylabel']
        if not options['yunit'] or options['yunit'] == options['heatmap.yunit']:
            options['yunit'] = options['diff.yunit']


        options['hide_y_ticklabels'] = True
        options['hide_y_ticks'] = True

    
        if options['plot_diff'] and len(data['path']) == 2:
            diff = data['diffractogram'][0]
            diff['I'] = diff['I'] - data['diffractogram'][1]['I']
            diff['I'] = diff['I'] - 0.5

            diff.plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=next(colours))


    # Adjust the plot to make it prettier
    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)
 

    # PLOT REFLECTION TABLES
    if options['reflections_plot'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        options['to_wavelength'] = data['wavelength'][0] # By default, the wavelength of the first diffractogram will be used for these.
        
        # Plot each reflection table in the relevant axis
        for reflections_params, axis in zip(options['reflections_data'], ref_axes):
            plot_reflection_table(data=data, reflections_params=reflections_params, ax=axis, options=options)

    # Print the reflection indices. 
    if options['reflections_indices'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        options['to_wavelength'] = data['wavelength'][0] # By default, the wavelength of the first diffractogram will be used for this.

        for reflections_params in options['reflections_data']:
            plot_reflection_indices(data=data, reflections_params=reflections_params, ax=indices_ax, options=options)


    ############################################################################################################################################################
    ##### UPDATE WIDGET ########################################################################################################################################
    ############################################################################################################################################################  
    
    if options['interactive_session_active']:
        options['current_y_offset'] = options['widget'].kwargs['offset_y']
        update_widgets(data=data, options=options)



    return data['diffractogram'], fig, ax



def generate_heatmap(data, options={}):

    required_options = ['x_tick_locators', 'heatmap_y_tick_locators', 'heatmap_normalise', 'normalisation_range', 'increase_contrast']

    default_options = {
        'x_tick_locators': [0.5, 0.1],
        'heatmap_y_tick_locators': [10, 5], # Major ticks for every 10 scans, minor for every 5
        'heatmap_normalise': False,
        'normalisation_range': None,
        'increase_contrast': False,
        'contrast_factor': 100
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    twotheta = []
    intensities = []
    scans = []

    for i, d in enumerate(data['diffractogram']):

        # Find normalisation factor
        if options['heatmap_normalise'] and options['normalisation_range']:
            mean_background = d['I'].loc[(d['2th'] > options['normalisation_range'][0]) & (d['2th'] < options['normalisation_range'][1])].mean()

            d['I'] = d['I'] / mean_background


        if options['increase_contrast']:

            if d['I'].min() < 0:
                d['I'] = d['I'] - d['I'].min()

            d['I'] = d['I']**(1/options['contrast_factor'])

        twotheta = np.append(twotheta, d['2th'].to_numpy())
        intensities = np.append(intensities, d['I'].to_numpy())
        scans = np.append(scans, np.full(len(d['2th'].to_numpy()), int(i)))
                

    heatmap = pd.DataFrame({'2th': twotheta, 'scan': scans, 'I': intensities})
    xrd.io.translate_wavelengths(data=heatmap, wavelength=data['wavelength'][0])

    min_dict = {'2th': heatmap['2th'].min(), '2th_cuka': heatmap['2th_cuka'].min(), '2th_moka': heatmap['2th_moka'].min(),
                        'q': heatmap['q'].min(), 'q2': heatmap['q2'].min(), 'q4': heatmap['q4'].min(), '1/d': heatmap['1/d'].min()}

    max_dict = {'2th': heatmap['2th'].max(), '2th_cuka': heatmap['2th_cuka'].max(), '2th_moka': heatmap['2th_moka'].max(),
                        'q': heatmap['q'].max(), 'q2': heatmap['q2'].max(), 'q4': heatmap['q4'].max(), '1/d': heatmap['1/d'].max()}


    ndatapoints = len(data['diffractogram'][0]['2th'])

    xlims = [0, ndatapoints, 0, ndatapoints] # 0: xmin, 1: xmax, 2: xmin_start, 3: xmax_start
    xticks = {}
    xticklabels = {}

    for xval in min_dict.keys():
       
       # Add xticks labels
        label_max = aux.floor(max_dict[xval], roundto=options['x_tick_locators'][0])
        label_min = aux.ceil(min_dict[xval], roundto=options['x_tick_locators'][0])
        label_steps = (label_max - label_min)/options['x_tick_locators'][0]

        xticklabels[xval] = np.linspace(label_min, label_max, num=int(label_steps)+1)

        # Add xticks
        xval_span = max_dict[xval] - min_dict[xval]
        steps = xval_span / ndatapoints
        
        
        xticks_xval = []

        for tick in xticklabels[xval]:
            xticks_xval.append((tick-min_dict[xval])/steps)

        xticks[xval] = xticks_xval


    options['x_tick_locators'] = None

    heatmap = heatmap.reset_index().pivot_table(index='scan', columns='2th', values='I')

    options['heatmap_xlim'] = xlims



    # Get temperatures if HTXRD-scans
    scan_numbers = []
    
    temperatures = []
    
    # FIXME This is a very bad check for whether it is HTXRD or not - it bascailly just excludes any files that has a .poni-file passed. Make more rigorous in future!
    if not 'calibrant' in data.keys():
        for i, filename in enumerate(data['path']):
            scan_numbers.append(i)
            temperatures.append(xrd.io.read_metadata_from_xy(filename)['temperature'])

        yticks = scan_numbers[0::options['heatmap_y_tick_locators'][0]]
        yticks.append(scan_numbers[-1])

        if not temperatures[0]:
            yticklabels = yticks
            options['heatmap.ylabel'] = 'Scan number'
            options['heatmap.yunit'] = None

        else:
            yticklabels = temperatures[0::options['heatmap_y_tick_locators'][0]]
            yticklabels.append(temperatures[-1])
            options['heatmap.ylabel'] = 'Temperature'
            options['heatmap.yunit'] = '$^\circ$C'

    else:
        yticks, yticklabels = None, None




    return heatmap, xticks, xticklabels, yticks, yticklabels







# #results = np.transpose(np.vstack([twotheta, scans, intensities]))


def determine_grid_layout(options):


    #aspect_ratio = int(options['format_params']['aspect_ratio'].split(':')[0]) / int(options['format_params']['aspect_ratio'].split(':')[1])

    nrows = 1 if not options['reflections_indices'] else 2

    if options['reflections_plot']:
        for reference in options['reflections_data']:
            nrows += 1

    options['format_params']['nrows'] = nrows

    if not 'grid_ratio_height' in options['format_params'].keys():
        options['format_params']['grid_ratio_height'] = [0.6 for i in range(nrows-1)]+[10]

    return options





def plot_diffractogram_interactive(data, options):


    # Format here is xminmax[0]: xmin, xminmax[1]: xmax, xminmax[2]: xmin_start, xminmax[3]: xmax_start, where "_start" denotes starting value of the slider
    xminmax = { '2th':      [None, None, None, None], '2th_cuka':   [None, None, None, None],   '2th_moka':     [None, None, None, None], 
                'd':        [None, None, None, None], '1/d':        [None, None, None, None], 
                'q':        [None, None, None, None], 'q2':         [None, None, None, None],   'q4':           [None, None, None, None], 
                'heatmap':  [None, None, None, None], 'start':      [None, None, None, None]}

    yminmax = { 'diff':     [None, None, None, None], 'heatmap':    [None, None, None, None],   'start':        [None, None, None, None]}
    
    update_xminmax(xminmax=xminmax, data=data, options=options)
    update_yminmax(yminmax=yminmax, data=data, options=options)

    options['xminmax'], options['yminmax'] = xminmax, yminmax

    # Get start values for ylim slider based on choice (FIXME This can be impleneted into update_yminmax). Can also make a 'start' item that stores the start values, instead of having 4 items in 'diff' as it is now.
    if options['heatmap']:
        ymin = yminmax['heatmap'][0]
        ymax = yminmax['heatmap'][1]
        ymin_start = yminmax['heatmap'][0]
        ymax_start = yminmax['heatmap'][1]

    elif not options['heatmap']:
        ymin = yminmax['diff'][0]
        ymax = yminmax['diff'][1]
        ymin_start = yminmax['diff'][2]
        ymax_start = yminmax['diff'][3]

    
    # FIXME The start values for xlim should probably also be decided by initial value of x_vals, and can likewise be implemented in update_xminmax() 


        
    options['widgets'] = {
        'xlim': {
            'w': widgets.FloatRangeSlider(value=[xminmax['start'][2], xminmax['start'][3]], min=xminmax['start'][0], max=xminmax['start'][1], step=0.5, layout=widgets.Layout(width='95%')),
            'state': options['x_vals'],
            '2th_default':      {'min': xminmax['2th'][0],          'max': xminmax['2th'][1],           'value': [xminmax['2th'][0],            xminmax['2th'][1]],             'step': 0.5},
            '2th_cuka_default': {'min': xminmax['2th_cuka'][0],     'max': xminmax['2th_cuka'][1],      'value': [xminmax['2th_cuka'][0],       xminmax['2th_cuka'][1]],        'step': 0.5},
            '2th_moka_default': {'min': xminmax['2th_moka'][0],     'max': xminmax['2th_moka'][1],      'value': [xminmax['2th_moka'][0],       xminmax['2th_moka'][1]],        'step': 0.5},
            'd_default':        {'min': xminmax['d'][0],            'max': xminmax['d'][1],             'value': [xminmax['d'][0],              xminmax['d'][1]],               'step': 0.5},
            '1/d_default':      {'min': xminmax['1/d'][0],          'max': xminmax['1/d'][1],           'value': [xminmax['1/d'][0],            xminmax['1/d'][1]],             'step': 0.5},
            'q_default':        {'min': xminmax['q'][0],            'max': xminmax['q'][1],             'value': [xminmax['q'][0],              xminmax['q'][1]],               'step': 0.5},
            'q2_default':       {'min': xminmax['q2'][0],           'max': xminmax['q2'][1],            'value': [xminmax['q2'][0],             xminmax['q2'][1]],              'step': 0.5},
            'q4_default':       {'min': xminmax['q4'][0],           'max': xminmax['q4'][1],            'value': [xminmax['q4'][0],             xminmax['q4'][1]],              'step': 0.5},
            'heatmap_default':  {'min': xminmax['heatmap'][0],      'max': xminmax['heatmap'][1],       'value': [xminmax['heatmap'][0],        xminmax['heatmap'][1]],         'step': 10}
        },
        'ylim': { 
            'w': widgets.FloatRangeSlider(value=[yminmax['start'][2], yminmax['start'][3]], min=yminmax['start'][0], max=yminmax['start'][1], step=0.01, layout=widgets.Layout(width='95%')),
            'state': 'heatmap' if options['heatmap'] else 'diff',
            'diff_default':     {'min': yminmax['diff'][0],         'max': yminmax['diff'][1],          'value': [yminmax['diff'][2],           yminmax['diff'][3]],            'step': 0.01},
            'heatmap_default':  {'min': yminmax['heatmap'][0],      'max': yminmax['heatmap'][1],       'value': [yminmax['heatmap'][0],        yminmax['heatmap'][1]],         'step': 0.01}
        }
    }

    if options['reflections_data']:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), data=widgets.fixed(data), options=widgets.fixed(options), 
        x_vals=widgets.Dropdown(options=['2th', 'd', '1/d', 'q', 'q2', 'q4', '2th_cuka', '2th_moka'], value='2th', description='X-values'),
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True), 
        reflections_plot=widgets.ToggleButton(value=True),
        reflections_indices=widgets.ToggleButton(value=False),
        heatmap=widgets.ToggleButton(value=options['heatmap']),
        xlim=options['widgets']['xlim']['w'],
        ylim=options['widgets']['ylim']['w'],
        offset_y=widgets.BoundedFloatText(value=options['offset_y'], min=-5, max=5, step=0.01, description='offset_y'),
        offset_x=widgets.BoundedFloatText(value=options['offset_x'], min=-1, max=1, step=0.01, description='offset_x')
        )
    
    else:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), data=widgets.fixed(data), options=widgets.fixed(options), 
        x_vals=widgets.Dropdown(options=['2th', 'd', '1/d', 'q', 'q2', 'q4', '2th_cuka', '2th_moka'], value='2th', description='X-values'),
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True),
        heatmap=widgets.ToggleButton(value=options['heatmap']),
        xlim=options['widgets']['xlim']['w'],
        ylim=options['widgets']['ylim']['w'],
        offset_y=widgets.BoundedFloatText(value=options['offset_y'], min=-5, max=5, step=0.01, description='offset_y'),
        offset_x=widgets.BoundedFloatText(value=options['offset_x'], min=-1, max=1, step=0.01, description='offset_x'))
    
    
    options['widget'] = w

    display(w)


def update_xminmax(xminmax, data, options={}):
    ''' Finds minimum and maximum values of each column and updates the minmax dictionary to contain the correct values.
    
    Input:
    minmax (dict): contains '''

    xminmax['2th'] = [None, None, None, None]
    for index, diffractogram in enumerate(data['diffractogram']):
        
        if not xminmax['2th'][0] or diffractogram['2th'].min() < xminmax['2th'][0]:
            xminmax['2th'][0] = diffractogram['2th'].min()
            min_index = index

        if not xminmax['2th'][1] or diffractogram['2th'].max() > xminmax['2th'][1]:
            xminmax['2th'][1] = diffractogram['2th'].max()
            max_index = index


    xminmax['2th'][2], xminmax['2th'][3] = xminmax['2th'][0], xminmax['2th'][1]

    xminmax['2th_cuka'][0], xminmax['2th_cuka'][1]  = data['diffractogram'][min_index]['2th_cuka'].min(),  data['diffractogram'][max_index]['2th_cuka'].max()
    xminmax['2th_cuka'][2], xminmax['2th_cuka'][3]  = xminmax['2th_cuka'][0], xminmax['2th_cuka'][1]

    xminmax['2th_moka'][0], xminmax['2th_moka'][1]  = data['diffractogram'][min_index]['2th_moka'].min(),  data['diffractogram'][max_index]['2th_moka'].max()
    xminmax['2th_moka'][2], xminmax['2th_moka'][3]  = xminmax['2th_moka'][0], xminmax['2th_moka'][1]

    xminmax['d'][0], xminmax['d'][1]                = data['diffractogram'][max_index]['d'].min(),         data['diffractogram'][min_index]['d'].max() # swapped, intended
    xminmax['d'][2], xminmax['d'][3]                = xminmax['d'][0], xminmax['d'][1]

    xminmax['1/d'][0], xminmax['1/d'][1]            = data['diffractogram'][min_index]['1/d'].min(),       data['diffractogram'][max_index]['1/d'].max()
    xminmax['1/d'][2], xminmax['1/d'][3]            = xminmax['1/d'][0], xminmax['1/d'][1]

    xminmax['q'][0], xminmax['q'][1]                = data['diffractogram'][min_index]['q'].min(),         data['diffractogram'][max_index]['q'].max()
    xminmax['q'][2], xminmax['q'][3]                = xminmax['q'][0], xminmax['q'][1]

    xminmax['q2'][0], xminmax['q2'][1]              = data['diffractogram'][min_index]['q2'].min(),        data['diffractogram'][max_index]['q2'].max()
    xminmax['q2'][2], xminmax['q2'][3]              = xminmax['q2'][0], xminmax['q2'][1]

    xminmax['q4'][0], xminmax['q4'][1]              = data['diffractogram'][min_index]['q4'].min(),        data['diffractogram'][max_index]['q4'].max()
    xminmax['q4'][2], xminmax['q4'][3]              = xminmax['q4'][0], xminmax['q4'][1]


    xminmax['heatmap'] = options['heatmap_xlim'] # This value is set in the generate_heatmap()-function


    xminmax['start'][0], xminmax['start'][1] = xminmax[options['x_vals']][0], xminmax[options['x_vals']][1]
    xminmax['start'][2], xminmax['start'][3] = xminmax[options['x_vals']][2], xminmax[options['x_vals']][3]


def update_yminmax(yminmax: dict, data: dict, options={}) -> None:
    
    yminmax['diff'] = [None, None, None, None]
    # Go through diffractograms and find the minimum and maximum intensity values
    for diffractogram in data['diffractogram']:
        if not yminmax['diff'][0] or (yminmax['diff'][0] > (diffractogram['I'].min())): 
            yminmax['diff'][0] = diffractogram['I'].min()

        if not yminmax['diff'][1] or (yminmax['diff'][1] < (diffractogram['I'].max())):
            yminmax['diff'][1] = diffractogram['I'].max()


    # Set start values of ymin and ymax to be slightly below lowest data points and slightly above highest data points to give some whitespace around the plot
    yminmax['diff'][2] = yminmax['diff'][0] - 0.1*yminmax['diff'][1]
    yminmax['diff'][3] = yminmax['diff'][1] + 0.2*yminmax['diff'][1]
    
    # Allow for adjustment up to five times ymax above and below data
    yminmax['diff'][0] = yminmax['diff'][0] - 5*yminmax['diff'][1]
    yminmax['diff'][1] = yminmax['diff'][1]*5


    # Set start values to the edges of the dataset
    yminmax['heatmap'][0], yminmax['heatmap'][1] = 0, data['heatmap'].shape[0]
    yminmax['heatmap'][2], yminmax['heatmap'][3] = yminmax['heatmap'][0], yminmax['heatmap'][1]


    if options['heatmap']:
        yminmax['start'][0], yminmax['start'][1] = yminmax['heatmap'][0], yminmax['heatmap'][1]
        yminmax['start'][2], yminmax['start'][3] = yminmax['heatmap'][0], yminmax['heatmap'][1]

    else:
        # The third and fourth index are different here to not be zoomed completely out to begin with.
        yminmax['start'][0], yminmax['start'][1] = yminmax['diff'][0], yminmax['diff'][1]
        yminmax['start'][2], yminmax['start'][3] = yminmax['diff'][2], yminmax['diff'][3]
    

def update_defaults(widget: dict, minmax: dict) -> None:
    ''' Updates the default x- or y-limits of a given widget. Refer to plot_diffractogram_interactive() to see the form of the widget that is passed in. An update of the min/max-values is done just prior to calling this function.
    Changes dictionaries in place.
    
    Input:
    widget (dict): A dictionary containing the widget itself (widget['w']) and all its default-values (e.g. widget['2th_default'])
    minmax (dict): A dictionary containing min and max values, as well as min_start and max_start values. (e.g. minmax['2th'] is a list with four elements: [xmin, xmax, xmin_start, xmax_start])
    
    Output:
    None.'''

    for name, attr in widget.items():
        if name.endswith('default'):
            attr['min'] = minmax[name.replace('_default', '')][0]
            attr['max'] = minmax[name.replace('_default', '')][1]           
            attr['value'] = [minmax[name.replace('_default', '')][2], minmax[name.replace('_default', '')][3]]


def update_widgets(data, options):


    for widget_name, widget in options['widgets'].items():

        # Make changes to xlim-widget
        if widget_name == 'xlim':
            # First update the min and max values
            update_xminmax(xminmax=options['xminmax'], data=data, options=options)
            update_defaults(widget=widget, minmax=options['xminmax'])
           

            if options['heatmap'] and (widget['state'] != 'heatmap'):


                setattr(widget['w'], 'min', widget['heatmap_default']['min'])
                setattr(widget['w'], 'max', widget['heatmap_default']['max'])
                setattr(widget['w'], 'value', widget['heatmap_default']['value'])
                setattr(widget['w'], 'step', widget['heatmap_default']['step'])

                widget['state'] = 'heatmap'
            
            elif not options['heatmap'] and (widget['state'] != options['x_vals']):              
                # Then loop through all attributes in the widget and change to current mode.
                for arg in widget[f'{options["x_vals"]}_default']:
                    
                    # If new min value is larger than previous max, or new max value is smaller than previous min, set the opposite first
                    if arg == 'min':
                        if widget[f'{options["x_vals"]}_default']['min'] > getattr(widget['w'], 'max'):
                            setattr(widget['w'], 'max', widget[f'{options["x_vals"]}_default']['max'])
                    
                    elif arg == 'max':
                        if widget[f'{options["x_vals"]}_default']['max'] < getattr(widget['w'], 'min'):
                            setattr(widget['w'], 'min', widget[f'{options["x_vals"]}_default']['min'])

                
                    setattr(widget['w'], arg, widget[f'{options["x_vals"]}_default'][arg])
                
                
                widget['state'] = options['x_vals']

        # Make changes to ylim-widget
        elif widget_name == 'ylim':
            update_yminmax(yminmax=options['yminmax'], data=data, options=options)
            update_defaults(widget=widget, minmax=options['yminmax'])
            
            state = 'heatmap' if options['heatmap'] else 'diff' 
    
            if widget['state'] != state or options['offset_change']:

                for arg in widget[f'{state}_default']:
                    # If new min value is larger than previous max, or new max value is smaller than previous min, set the opposite first
                    if arg == 'min':                     
                        if widget[f'{state}_default']['min'] > getattr(widget['w'], 'max'):
                            setattr(widget['w'], 'max', widget[f'{state}_default']['max'])
                    
                    elif arg == 'max':
                        if widget[f'{state}_default']['max'] < getattr(widget['w'], 'min'):
                            setattr(widget['w'], 'min', widget[f'{state}_default']['min'])

                
                    setattr(widget['w'], arg, widget[f'{state}_default'][arg])
                
                options['offset_change'] = False
                widget['state'] = state




def plot_reflection_indices(data, reflections_params, ax, options={}):
    ''' Print reflection indices from output generated by VESTA.
    
    Required contents of data:
    path (str): relative path to reflection table file'''

    required_options = ['reflection_indices', 'text_colour', 'hide_indices']

    default_options = {
        'reflection_indices': 3, # Number of reflection indices to plot, from highest intensity and working its way down
        'text_colour': 'black',
        'hide_indices': False
    }

    reflections_params = aux.update_options(options=reflections_params, required_options=required_options, default_options=default_options)

    if not reflections_params['hide_indices']:
        reflection_table = xrd.io.load_reflection_table(data=data, reflections_params=reflections_params, options=options)
        
        if reflections_params['reflection_indices'] > 0:

            # Get the data['reflection_indices'] number of highest reflections within the subrange options['xlim']
            x_vals = 'heatmap' if options['heatmap'] else options['x_vals']
            reflection_indices = reflection_table.loc[(reflection_table[x_vals] > options['xlim'][0]) & (reflection_table[x_vals] < options['xlim'][1])].nlargest(options['reflection_indices'], 'I')

            # Plot the indices
            for i in range(reflections_params['reflection_indices']):
                if reflection_indices.shape[0] > i:
                    ax.text(s=f'({reflection_indices["h"].iloc[i]} {reflection_indices["k"].iloc[i]} {reflection_indices["l"].iloc[i]})', x=reflection_indices[x_vals].iloc[i], y=0, fontsize=2.5, rotation=90, va='bottom', ha='center', c=reflections_params['text_colour'])    

        
        if options['xlim']:
            ax.set_xlim(options['xlim'])
        
        ax.axis('off')


    return
    
def plot_reflection_table(data, reflections_params, ax=None, options={}):
    ''' Plots a reflection table from output generated by VESTA.
    
    Required contents of data:
    path (str): relative path to reflection table file'''

    required_options = ['reflection_indices', 'reflections_colour', 'min_alpha', 'wavelength', 'format_params', 'rc_params', 'label']

    default_options = {
        'reflection_indices': 0, # Number of indices to print
        'reflections_colour': [0,0,0],
        'min_alpha': 0,
        'wavelength': 1.54059, # CuKalpha, [Å] 
        'format_params': {},
        'rc_params': {},
        'label': None
    }


    if 'colour' in reflections_params.keys():
        options['reflections_colour'] = reflections_params['colour']
    if 'min_alpha' in reflections_params.keys():
        options['min_alpha'] = reflections_params['min_alpha']
    if 'reflection_indices' in reflections_params.keys():
        options['reflection_indices'] = reflections_params['reflection_indices']
    if 'label' in reflections_params.keys():
        options['label'] = reflections_params['label']
    if 'wavelength' in reflections_params.keys():
        options['wavelength'] = reflections_params['wavelength']

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    
    if not ax:
        _, ax = btp.prepare_plot(options)

    x_vals = 'heatmap' if options['heatmap'] else options['x_vals']

    reflection_table = xrd.io.load_reflection_table(data=data, reflections_params=reflections_params, options=options)
    reflections, intensities  = reflection_table[x_vals], reflection_table['I']


    
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

        ax.text(s=reflections_params['label'], x=(ax.get_xlim()[0]-0.01*xlim_range), y=ylim_avg, ha = 'right', va = 'center')





def prettify_labels(label):
	
	labels_dict = {
		'2th': '2$\\theta$',
        'I': 'Intensity'	
        }

	return labels_dict[label]



def reverse_diffractograms(diffractograms):

    rev_diffractograms = []

    for i in len(diffractograms):
        rev_diffractograms.append(diffractograms.pop())

    return rev_diffractograms



def make_animation(data: dict, options={}):

    default_options = {
        'cmap': 'inferno',
        'contrast': False,
        'contrast_factor': 1/3,
        'save_path': 'diff_animation.gif',
        'fps': 5
    }

    options = aux.update_options(options=options, default_options=default_options, required_options=default_options.keys())

    if not isinstance(data['path'], list):
        data['path'] = aux.get_filenames(data['path'], ext='dat')

    
    if not os.path.isdir('tmp'):
        os.makedirs('tmp')
    
	# Scale image to make GIF smaller
    # 
    options['format_params']['width'] = 5
    options['format_params']['height'] = 5

    options['format_params']['dpi'] = 200

    for i, scan in enumerate(data['path']):
            
        giffig, gifax = btp.prepare_plot(options=options)

        img = xrd.io.get_image_array(scan)
        
        if options['contrast']:
            img[img < 0] = 0.00000001
            img = np.log(img)
            img[img < 0] = 0

        gifax.imshow(img, cmap=options['cmap'])

        plt.savefig(os.path.join('tmp', str(i+1).zfill(4)+'.png'))
        plt.close()

        
    img_paths = aux.get_filenames('tmp', ext='png')
    
    frames = []
    for path in img_paths:
        frame = Image.open(path)
        frames.append(frame)

    frames[0].save(options['save_path'], format='GIF', append_images=frames[1:], save_all=True, duration=(1/options['fps'])*1000, loop=0)

    shutil.rmtree('tmp')



def plot_refinement(data, options={}):


    required_options = ['diff_offset', 'index', 'title', 'xlim', 'r_wp', 'r_exp', 'wp']

    default_options = {
        'diff_offset': .10,
        'index': -1,
        'title': None,
        'xlim': None,
        'r_wp': True,
        'r_exp': False,
        'wp': False,
    }

    options = aux.update_options(options=options, default_options=default_options, required_options=required_options)

    df = pd.read_csv(data['path'], delim_whitespace=True, header=None)
    df.columns = ['2th', 'Yobs', 'Ycalc', 'diff']
    df['diff'] = df['diff'] - options['diff_offset']*(df['Yobs'].max() - df['Yobs'].min())
    

    if not isinstance(data['results'], list):
        data['results'] = [data['results']]

    results = {
        'vol': [],
        'mass': [],
        'wp': [],
        'a': [],
        'b': [],
        'c': [],
        'alpha': [],
        'beta': [],
        'gamma': []
    }

    for result in data['results']:
        result = xrd.refinement.read_results(path=result)

        r_wp = result['r_wp'].iloc[options['index']]
        r_exp = result['r_exp'].iloc[options['index']]

        for attr in results.keys():
            results[attr].append(result[attr].iloc[options['index']])

    fig, ax = plt.subplots(figsize=(20,10))

    df.plot(x='2th', y='Yobs', kind='scatter', ax=ax, c='black', marker='$\u25EF$')
    df.plot(x='2th', y='Ycalc', ax=ax, c='red')
    df.plot(x='2th', y='diff', ax=ax)
    
    if options['r_wp']:
        ax.text(x=0.7*df['2th'].max(), y=0.7*df['Yobs'].max(), s='R$_{wp}$ = '+f'{r_wp}', fontsize=20)
    
    if options['r_exp']:
        ax.text(x=0.70*df['2th'].max(), y=0.60*df['Yobs'].max(), s='R$_{exp}$ = '+f'{r_exp}', fontsize=20)


    if options['wp']:
        for i, (result, label) in enumerate(zip(data['results'], options['labels'])):
            ax.text(x=0.7*df['2th'].max(), y=(0.9-0.1*i)*df['Yobs'].max(), s=f'{label}: {np.round(float(results["wp"][i]), 2)}%', fontsize=20)

    if options['title']:
        ax.set_title(options['title'], size=30)

    if options['xlim']:
        ax.set_xlim(options['xlim'])
    else:
        ax.set_xlim([df['2th'].min(), df['2th'].max()])

    ax.tick_params(which='both', labelleft=False, left=False, labelsize=20, direction='in')
    ax.set_ylabel('Intensity [arb. u.]', size=20)
    ax.set_xlabel('2$\\theta$ [$^{\circ}$]', size=20)