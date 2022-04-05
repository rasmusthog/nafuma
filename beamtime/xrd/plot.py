import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math

import ipywidgets as widgets

import beamtime.xrd as xrd
import beamtime.auxillary as aux
import beamtime.plotting as btp


def plot_diffractogram(data, options={}):
    ''' Plots a diffractogram.
    
    Input:
    data (dict): Must include path = string to diffractogram data, and plot_kind = (recx, beamline, image)'''

    # Update options
    required_options = ['x_vals', 'y_vals', 'ylabel', 'xlabel', 'xunit', 'yunit', 'line', 'scatter', 'xlim', 'ylim', 'normalise', 'offset', 'offset_x', 'offset_y', 'offset_change',
    'reflections_plot', 'reflections_indices', 'reflections_data', 'heatmap', 'cmap', 'plot_kind', 'palettes', 'interactive', 'rc_params', 'format_params', 'interactive_session_active']

    default_options = {
        'x_vals': '2th', 
        'y_vals': 'I',
        'ylabel': 'Intensity', 'xlabel': '2theta', 
        'xunit': 'deg', 'yunit': 'a.u.',
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
        'interactive': False,
        'interactive_session_active': False,
        'rc_params': {},
        'format_params': {},
        }

    if 'offset_y' not in options.keys():
        if len(data['path']) > 10:
            default_options['offset_y'] = 0.05

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    options['current_offset_y'] = options['offset_y']

    # Convert data['path'] to list to allow iteration over this to accommodate both single and multiple diffractograms
    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

   

    # Check if there is some data stored already, load in data if not. This speeds up replotting in interactive mode.
    if not 'diffractogram' in data.keys():
        # Initialise empty list for diffractograms and wavelengths
        data['diffractogram'] = [None for _ in range(len(data['path']))]
        data['wavelength'] = [None for _ in range(len(data['path']))]

        for index in range(len(data['path'])):
            diffractogram, wavelength = xrd.io.read_data(data=data, options=options, index=index)
            
            data['diffractogram'][index] = diffractogram
            data['wavelength'][index] = wavelength

        # Sets the xlim if this has not bee specified
        if not options['xlim']:
            options['xlim'] = [data['diffractogram'][0][options['x_vals']].min(), data['diffractogram'][0][options['x_vals']].max()]

        # Generate heatmap data
        data['heatmap'], data['heatmap_xticks'], data['heatmap_xticklabels'] = generate_heatmap(data=data, options=options)
        if options['heatmap']:
            options['xlim'] = options['heatmap_xlim']

    else:
        if not isinstance(data['diffractogram'], list):
            data['diffractogram'] = [data['diffractogram']]
            data['wavelength'] = [data['wavelength']]



    if options['interactive_session_active']:
        if options['offset']:
            if (options['offset_x'] != options['current_offset_x']) or (options['offset_y'] != options['current_offset_y']):
                for i, (diff, wl) in enumerate(zip(data['diffractogram'], data['wavelength'])):
                    xrd.io.apply_offset(diff, wl, i, options)


    # Start inteactive session with ipywidgets. Disables options['interactive'] in order for the interactive loop to not start another interactive session
    if options['interactive']:
        options['interactive'] = False
        options['interactive_session_active'] = True
        plot_diffractogram_interactive(data=data, options=options)
        return
    
    
    # Makes a list out of reflections_data if it only passed as a dict, as it will be looped through later
    if options['reflections_data']:
        if not isinstance(options['reflections_data'], list):
            options['reflections_data'] = [options['reflections_data']]

    # Determine number of subplots and height ratios between them
    if len(options['reflections_data']) >= 1:
        options = determine_grid_layout(options=options)


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

    if len(data['path']) < 10:
        colours = btp.generate_colours(options['palettes'])
    else:
        colours = btp.generate_colours(['black'], kind='single')

    if options['heatmap']:
        sns.heatmap(data['heatmap'], cmap=options['cmap'], cbar=False, ax=ax)
        ax.set_xticks(data['heatmap_xticks'][options['x_vals']])
        ax.set_xticklabels(data['heatmap_xticklabels'][options['x_vals']])
        ax.tick_params(axis='x', which='minor', bottom=False, top=False)

    else:
        for diffractogram in data['diffractogram']:
            if options['line']:
                diffractogram.plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=next(colours), zorder=1)
        
            if options['scatter']:
                ax.scatter(x=diffractogram[options['x_vals']], y = diffractogram[options['y_vals']], c=[(1,1,1,0)], edgecolors=[next(colours)], linewidths=plt.rcParams['lines.markeredgewidth'], zorder=2) #, edgecolors=np.array([next(colours)]))



    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)

    

    # Make the reflection plots. By default, the wavelength of the first diffractogram will be used for these.
    if options['reflections_plot'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        options['to_wavelength'] = data['wavelength'][0]
        
        for reference, axis in zip(options['reflections_data'], ref_axes):
            plot_reflection_table(data=reference, ax=axis, options=options)

    # Print the reflection indices. By default, the wavelength of the first diffractogram will be used for this.
    if options['reflections_indices'] and options['reflections_data']:
        options['xlim'] = ax.get_xlim()
        options['to_wavelength'] = data['wavelength'][0]

        for reference in options['reflections_data']:
            plot_reflection_indices(data=reference, ax=indices_ax, options=options)


    if options['interactive_session_active']:
        options['current_y_offset'] = options['widget'].kwargs['offset_y']
        update_widgets(data=data, options=options)



    return data['diffractogram'], fig, ax



def generate_heatmap(data, options={}):

    required_options = ['x_tick_locators']

    default_options = {
        'x_tick_locators': [0.5, 0.1]
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    twotheta = []
    intensities = []
    scans = []

    for i, d in enumerate(data['diffractogram']):
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


    return heatmap, xticks, xticklabels







# #results = np.transpose(np.vstack([twotheta, scans, intensities]))


def determine_grid_layout(options):


    nrows = 1 if not options['reflections_indices'] else 2

    if options['reflections_plot']:
        for reference in options['reflections_data']:
            nrows += 1

    options['format_params']['nrows'] = nrows
    options['format_params']['grid_ratio_height'] = [1 for i in range(nrows-1)]+[10]

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
            'w': widgets.FloatRangeSlider(value=[yminmax['start'][2], yminmax['start'][3]], min=yminmax['start'][0], max=yminmax['start'][1], step=0.5, layout=widgets.Layout(width='95%')),
            'state': 'heatmap' if options['heatmap'] else 'diff',
            'diff_default':     {'min': yminmax['diff'][0],         'max': yminmax['diff'][1],          'value': [yminmax['diff'][2],           yminmax['diff'][3]],            'step': 0.1},
            'heatmap_default':  {'min': yminmax['heatmap'][0],      'max': yminmax['heatmap'][1],       'value': [yminmax['heatmap'][0],        yminmax['heatmap'][1]],         'step': 0.1}
        }
    }

    if options['reflections_data']:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), data=widgets.fixed(data), options=widgets.fixed(options), 
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True), 
        reflections_plot=widgets.ToggleButton(value=True),
        reflections_indices=widgets.ToggleButton(value=False),
        heatmap=widgets.ToggleButton(value=options['heatmap']),
        x_vals=widgets.Dropdown(options=['2th', 'd', '1/d', 'q', 'q2', 'q4', '2th_cuka', '2th_moka'], value='2th', description='X-values'),
        xlim=options['widgets']['xlim']['w'],
        ylim=options['widgets']['ylim']['w'],
        offset_y=widgets.BoundedFloatText(value=options['offset_y'], min=-5, max=5, step=0.01, description='offset_y'),
        offset_x=widgets.BoundedFloatText(value=options['offset_x'], min=-1, max=1, step=0.01, description='offset_x')
        )
    
    else:
        w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_diffractogram), data=widgets.fixed(data), options=widgets.fixed(options), 
        scatter=widgets.ToggleButton(value=False), 
        line=widgets.ToggleButton(value=True),
        xlim=options['widgets']['xlim']['w'])
    
    
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




def plot_reflection_indices(data, ax, options={}):
    ''' Print reflection indices from output generated by VESTA.
    
    Required contents of data:
    path (str): relative path to reflection table file'''

    required_options = ['reflection_indices', 'text_colour', 'hide_indices']

    default_options = {
        'reflection_indices': 3, # Number of reflection indices to plot, from highest intensity and working its way down
        'text_colour': 'black',
        'hide_indices': False
    }

    data = aux.update_options(options=data, required_options=required_options, default_options=default_options)

    if not data['hide_indices']:
        reflection_table = xrd.io.load_reflection_table(data=data, options=options)
        
        if data['reflection_indices'] > 0:

            # Get the data['reflection_indices'] number of highest reflections within the subrange options['xlim']
            reflection_indices = reflection_table.loc[(reflection_table[options['x_vals']] > options['xlim'][0]) & (reflection_table[options['x_vals']] < options['xlim'][1])].nlargest(options['reflection_indices'], 'I')

            # Plot the indices
            for i in range(data['reflection_indices']):
                if reflection_indices.shape[0] > i:
                    ax.text(s=f'({reflection_indices["h"].iloc[i]} {reflection_indices["k"].iloc[i]} {reflection_indices["l"].iloc[i]})', x=reflection_indices[options['x_vals']].iloc[i], y=0, fontsize=2.5, rotation=90, va='bottom', ha='center', c=data['text_colour'])    

        
        if options['xlim']:
            ax.set_xlim(options['xlim'])
        
        ax.axis('off')


    return
    
def plot_reflection_table(data, ax=None, options={}):
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

    if 'colour' in data.keys():
        options['reflections_colour'] = data['colour']
    if 'min_alpha' in data.keys():
        options['min_alpha'] = data['min_alpha']
    if 'reflection_indices' in data.keys():
        options['reflection_indices'] = data['reflection_indices']
    if 'label' in data.keys():
        options['label'] = data['label']
    if 'wavelength' in data.keys():
        options['wavelength'] = data['wavelength']

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    
    if not ax:
        _, ax = btp.prepare_plot(options)

    reflection_table = xrd.io.load_reflection_table(data=data, options=options)

    reflections, intensities  = reflection_table[options['x_vals']], reflection_table['I']


    
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

        ax.text(s=data['label'], x=(ax.get_xlim()[0]-0.01*xlim_range), y=ylim_avg, ha = 'right', va = 'center')





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
