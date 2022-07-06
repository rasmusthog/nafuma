from logging import raiseExceptions
from jinja2 import TemplateRuntimeError
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import nafuma.auxillary as aux
import nafuma.plotting as btp
import nafuma.xanes as xas
import nafuma.xanes.io as io
from scipy.signal import savgol_filter
from datetime import datetime
import ipywidgets as widgets
from IPython.display import display


##Better to make a new function that loops through the files, and performing the split_xanes_scan on

#Trying to make a function that can decide which edge it is based on the first ZapEnergy-value
def find_element(data: dict, index=0) -> str:
    ''' Takes the data dictionary and determines based on the start value of the ZapEnergy-column which element the edge is from.'''

    element_energy_intervals = {
        'Mn': [5.9, 6.5],
        'Fe': [7.0, 7.2],
        'Co': [7.6, 7.8],
        'Ni': [8.0, 8.6]
    }

    if (element_energy_intervals['Mn'][0] < data['xanes_data_original']["ZapEnergy"].iloc[index]) & (data['xanes_data_original']["ZapEnergy"].iloc[index] < element_energy_intervals['Mn'][1]):
        edge = 'Mn'
    elif (element_energy_intervals['Fe'][0] < data['xanes_data_original']["ZapEnergy"].iloc[index]) & (data['xanes_data_original']["ZapEnergy"].iloc[index] < element_energy_intervals['Fe'][1]):
        edge = 'Fe'
    elif (element_energy_intervals['Co'][0] < data['xanes_data_original']["ZapEnergy"].iloc[index]) & (data['xanes_data_original']["ZapEnergy"].iloc[index] < element_energy_intervals['Co'][1]):
        edge = 'Co'   
    elif (element_energy_intervals['Ni'][0] < data['xanes_data_original']["ZapEnergy"].iloc[index]) & (data['xanes_data_original']["ZapEnergy"].iloc[index] < element_energy_intervals['Ni'][1]):
        edge = 'Ni'
        
        
    return(edge)



def pre_edge_fit(data: dict, options={}) -> pd.DataFrame:


    # FIXME Add log-file

    required_options = ['pre_edge_limits', 'pre_edge_masks', 'pre_edge_polyorder', 'pre_edge_store_data', 'log', 'logfile', 'show_plots', 'save_plots', 'save_folder', 'ylim', 'interactive', 'interactive_session_active']
    default_options = {
        'pre_edge_limits': [None, None],
        'pre_edge_masks': [],
        'pre_edge_polyorder': 1,
        'pre_edge_store_data': False,
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_pre_edge_fit.log',
        'show_plots': False,
        'save_plots': False,
        'save_folder': './',
        'ylim': [None, None],
        'interactive': False,
        'interactive_session_active': False
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if options['log']:
        aux.write_log(message='Starting pre edge fit', options=options)

    # FIXME Implement with finding accurate edge position
    # FIXME Allow specification of start of pre-edge area
    # Find the cutoff point at which the edge starts - everything to the LEFT of this point will be used in the pre edge function fit
    if not options['pre_edge_limits'][0]:
        options['pre_edge_limits'][0] = data['xanes_data_original']['ZapEnergy'].min()

    
    if not options['pre_edge_limits'][1]:
        pre_edge_limit_offset = 0.03

        data['edge'] = find_element(data)

        edge_position = estimate_edge_position(data, options, index=0)
        options['pre_edge_limits'][1] = edge_position - pre_edge_limit_offset

        print(edge_position)

    if options['pre_edge_limits'][0] >= options['pre_edge_limits'][1]:
        options['pre_edge_limits'][1] = options['pre_edge_limits'][0] + 0.03

    # Start inteactive session with ipywidgets. Disables options['interactive'] in order for the interactive loop to not start another interactive session
    if options['interactive']:
        options['interactive'] = False
        options['interactive_session_active'] = True
        options['show_plots'] = True
        pre_edge_fit_interactive(data=data, options=options)
        return



    # FIXME There should be an option to specify the interval in which to fit the background - now it is taking everything to the left of edge_start parameter, but if there are some artifacts in this area, it should be possible to
    # limit the interval
    # Making a dataframe only containing the rows that are included in the background subtraction (points lower than where the edge start is defined)
    pre_edge_data = data['xanes_data_original'].loc[(data['xanes_data_original']["ZapEnergy"] > options['pre_edge_limits'][0]) & (data['xanes_data_original']["ZapEnergy"] < options['pre_edge_limits'][1])].copy()

    for mask in options['pre_edge_masks']:
        pre_edge_data.loc[(pre_edge_data['ZapEnergy'] > mask[0]) & (pre_edge_data['ZapEnergy'] < mask[1])] = np.nan

    pre_edge_data = pre_edge_data.dropna()

    # Making a new dataframe, with only the ZapEnergies as the first column -> will be filled to include the background data
    pre_edge_fit_data = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    data['pre_edge_params'] = {}

    for i, filename in enumerate(data['path']):

        if options['interactive_session_active'] and i > 0:
            continue

        if options['log']:
            aux.write_log(message=f'... Fitting pre edge on {os.path.basename(filename)} ({i+1}/{len(data["path"])})', options=options)

        #Fitting linear function to the background
        params = np.polyfit(pre_edge_data["ZapEnergy"],pre_edge_data[filename],options['pre_edge_polyorder'])
        fit_function = np.poly1d(params)

        data['pre_edge_params'][filename] = params

        if options['log']:
            aux.write_log(message=f'...... Pre edge fitted between {options["pre_edge_limits"][0]} and {options["pre_edge_limits"][1]} with polynomial of order {options["pre_edge_polyorder"]} with parmameters {params}.', options=options)
            if options['pre_edge_masks']:
                aux.write_log(message=f'...... Excluded regions: {options["pre_edge_masks"]}', options=options)
        
        #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        background=fit_function(pre_edge_fit_data["ZapEnergy"])
            
        #adding a new column in df_background with the y-values of the background
        pre_edge_fit_data.insert(1,filename,background) 
        
        if options['show_plots'] or options['save_plots']:
            fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,10))
            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax1)
            pre_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
            ax1.axvline(x = max(pre_edge_data['ZapEnergy']), ls='--')
            ax1.axvline(x = min(pre_edge_data['ZapEnergy']), ls='--')
            ax1.set_title(f'{os.path.basename(filename)} - Full view', size=20)
            
            if options['ylim'][0] != None:
                ax1.set_ylim(bottom=options['ylim'][0])
            if options['ylim'][1]:
                ax1.set_ylim(top=options['ylim'][1])

            for mask in options['pre_edge_masks']:
                ax1.fill_between(x=mask, y1=0, y2=data['xanes_data_original'][filename].max()*2, alpha=0.2, color='black')

            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax2)
            pre_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
            ax2.axvline(x = max(pre_edge_data['ZapEnergy']), ls='--')
            ax2.set_xlim([min(pre_edge_data['ZapEnergy']), max(pre_edge_data['ZapEnergy'])])
            ax2.set_ylim([min(pre_edge_data[filename]), max(pre_edge_data[filename])])
            ax2.set_title(f'{os.path.basename(filename)} - Fit region', size=20)

            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_pre_edge_fit.png'
                plt.savefig(dst, transparent=False)
            
            if not options['show_plots']:
                plt.close()


    if options['log']:
        aux.write_log(message=f'Pre edge fitting done.', options=options)

    if options['pre_edge_store_data']:
        data['pre_edge_fit_data'] = pre_edge_fit_data

    return pre_edge_fit_data



def pre_edge_fit_interactive(data: dict, options: dict) -> None:


    w = widgets.interactive(
        btp.ipywidgets_update, func=widgets.fixed(pre_edge_fit), data=widgets.fixed(data), options=widgets.fixed(options), 
        pre_edge_limits=widgets.FloatRangeSlider(value=[options['pre_edge_limits'][0], options['pre_edge_limits'][1]], min=data['xanes_data_original']['ZapEnergy'].min(), max=data['xanes_data_original']['ZapEnergy'].max(), step=0.0001),
        pre_edge_store_data=widgets.Checkbox(value=options['pre_edge_store_data'])
    )
    
    options['widget'] = w

    display(w)




def pre_edge_subtraction(data: dict, options={}):

    required_options = ['log', 'logfile', 'show_plots', 'save_plots', 'save_folder', 'pre_edge_subtraction_store_data']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_pre_edge_subtraction.log',
        'show_plots': False,
        'save_plots': False,
        'save_folder': './',
        'pre_edge_subtraction_store_data': True
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    if options['log']:
        aux.write_log(message='Starting pre edge subtraction', options=options)

    xanes_data_bkgd_subtracted = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    for i, filename in enumerate(data['path']):
        if options['log']:
            aux.write_log(message=f'... Subtracting background on {os.path.basename(filename)} ({i}/{len(data["path"])})', options=options)

        xanes_data_bkgd_subtracted.insert(1, filename, data['xanes_data_original'][filename] - data['pre_edge_fit_data'][filename])

        if options['save_plots'] or options['show_plots']:


                fig, ax = plt.subplots(figsize=(10,5))
                data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax, label='Original data')
                xanes_data_bkgd_subtracted.plot(x='ZapEnergy', y=filename, color='red', ax=ax, label='Pre edge subtracted')
                ax.set_title(f'{os.path.basename(filename)} - After subtraction', size=20)

                
                if options['save_plots']:
                    if not os.path.isdir(options['save_folder']):
                        os.makedirs(options['save_folder'])

                    dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_pre_edge_subtraction.png'

                    plt.savefig(dst)

                if not options['show_plots']:    
                    plt.close()

    if options['pre_edge_subtraction_store_data']:
        data['xanes_data'] = xanes_data_bkgd_subtracted

    return xanes_data_bkgd_subtracted





def post_edge_fit(data: dict, options={}):
    ''' Fit the post edge within the post_edge.limits to a polynomial of post_edge.polyorder order. Allows interactive plotting, as well as showing static plots and saving plots to drive.
    
    Requires data to have already been read to data['xanes_data_original'] 
    '''
    
    
    required_options = ['log', 'logfile', 'post_edge_masks', 'post_edge_limits', 'post_edge_polyorder', 'post_edge_store_data', 'interactive', 'interactive_session_active', 'show_plots', 'save_plots', 'save_folder']
    default_options = {
        'log': False, 
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_post_edge_fit.log',
        'post_edge_limits': [None, None],
        'post_edge_masks': [],
        'post_edge_polyorder': 2,
        'post_edge_store_data': False,
        'interactive': False,
        'interactive_session_active': False,
        'show_plots': False,
        'save_plots': False,
        'save_folder': './',
        'ylim': [None, None]
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if not options['post_edge_limits'][0]:
        post_edge_limit_offset = 0.03

        data['edge'] = find_element(data)

        edge_position = estimate_edge_position(data, options, index=0)
        options['post_edge_limits'][0] = edge_position + post_edge_limit_offset

    if not options['post_edge_limits'][1]:
        options['post_edge_limits'][1] = data['xanes_data_original']['ZapEnergy'].max()

    if options['post_edge_limits'][0] > options['post_edge_limits'][1]:
        options['post_edge_limits'][0] = options['post_edge_limits'][1] - 0.1
    

    # Start inteactive session with ipywidgets. Disables options['interactive'] in order for the interactive loop to not start another interactive session
    if options['interactive']:
        options['interactive'] = False
        options['interactive_session_active'] = True
        options['show_plots'] = True
        post_edge_fit_interactive(data=data, options=options)
        return



    post_edge_data = data['xanes_data_original'].loc[(data['xanes_data_original']["ZapEnergy"] > options['post_edge_limits'][0]) & (data['xanes_data_original']["ZapEnergy"] < options['post_edge_limits'][1])].copy()

    for mask in options['post_edge_masks']:
        post_edge_data.loc[(post_edge_data['ZapEnergy'] > mask[0]) & (post_edge_data['ZapEnergy'] < mask[1])] = np.nan

    post_edge_data = post_edge_data.dropna() #Removing all indexes without any value, as some of the data sets misses the few last data points and fucks up the fit



    # Making a new dataframe, with only the ZapEnergies as the first column -> will be filled to include the background data
    post_edge_fit_data = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    data['post_edge_params'] = {}
    
    for i, filename in enumerate(data['path']):

        if options['interactive_session_active'] and i > 0:
            continue

        if options['log']:
            aux.write_log(message=f'... Fitting post edge on {os.path.basename(filename)} ({i+1}/{len(data["path"])}) with polynomial order {options["post_edge_polyorder"]}', options=options)

        #Fitting linear function to the background
        params = np.polyfit(post_edge_data["ZapEnergy"], post_edge_data[filename], options['post_edge_polyorder'])
        fit_function = np.poly1d(params)

        if options['log']:
            aux.write_log(message=f'...... Post edge fitted between {options["post_edge_limits"][0]} and {options["post_edge_limits"][1]} with polynomial of order {options["post_edge_polyorder"]} with parmameters {params}.', options=options)
            if options['post_edge_masks']:
                aux.write_log(message=f'...... Excluded regions: {options["post_edge_masks"]}', options=options)

        data['post_edge_params'][filename] = params
        
        #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        background=fit_function(post_edge_fit_data["ZapEnergy"])
            
        #adding a new column in df_background with the y-values of the background
        post_edge_fit_data.insert(1,filename,background) 
        
        if options['save_plots'] or options['show_plots']:


            fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,10))
            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax1)
            post_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
            ax1.axvline(x = max(post_edge_data['ZapEnergy']), ls='--')
            ax1.axvline(x = min(post_edge_data['ZapEnergy']), ls='--')
            ax1.set_title(f'{os.path.basename(filename)} - Full view', size=20)

            for mask in options['post_edge_masks']:
                ax1.fill_between(x=mask, y1=0, y2=data['xanes_data_original'][filename].max()*2, alpha=0.2, color='black')

            if options['ylim'][0] != None:
                ax1.set_ylim(bottom=options['ylim'][0])
            if options['ylim'][1] != None:
                ax1.set_ylim(top=options['ylim'][1])

            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax2)
            post_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
            ax2.axvline(x = max(post_edge_data['ZapEnergy']), ls='--')
            ax2.set_xlim([min(post_edge_data['ZapEnergy']), max(post_edge_data['ZapEnergy'])])
            ax2.set_ylim([min(post_edge_data[filename]), max(post_edge_data[filename])])
            ax2.set_title(f'{os.path.basename(filename)} - Fit region', size=20)

            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_post_edge_fit.png'

                plt.savefig(dst, transparent=False)
            
            if not options['show_plots']:
                plt.close()


    if options['log']:
        aux.write_log(message='Post edge fitting done!', options=options)

    if options['post_edge_store_data']:
        data['post_edge_fit_data'] = post_edge_fit_data.dropna(axis=0)


    return post_edge_fit_data


def post_edge_fit_interactive(data: dict, options: dict) -> None:
    ''' Defines the widgets to use with the ipywidgets interactive mode and calls the update function found in btp.ipywidgets. '''

    w = widgets.interactive(
        btp.ipywidgets_update, func=widgets.fixed(post_edge_fit), data=widgets.fixed(data), options=widgets.fixed(options), 
        post_edge_limits=widgets.FloatRangeSlider(value=[options['post_edge_limits'][0], options['post_edge_limits'][1]], min=data['xanes_data_original']['ZapEnergy'].min(), max=data['xanes_data_original']['ZapEnergy'].max(), step=0.0001),
        post_edge_store_data=widgets.Checkbox(value=options['post_edge_store_data'])
    )
    
    options['widget'] = w

    display(w)

def smoothing(data: dict, options={}):
    ' Smoothes the data using the Savitzky-Golay filter. This is the only algorithm at this moment.  '


    required_options = ['log', 'logfile', 'show_plots', 'save_plots', 'save_folder', 'interactive', 'smooth_window_length', 'smooth_algorithm', 'smooth_polyorder', 'smooth_save_default', 'smooth_store_data']
    default_options = {
        'log': False, # Toggles logging on / off
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_smoothing.log', # Sets path to log-file. Ignored if log == False
        'show_plots': False, # Toggles showing plots on / off. This is only recommended when working with a handful of scans. 
        'save_plots': False, # Toggles saving plots on / off
        'save_folder': './', # Sets path to folder where plots should be saved. Ignored if save_plots == False
        'interactive': False, # Toggles interactive mode on / off. This is only recommended for a single scan to determine proper parameters for smoothing.
        'smooth_window_length': 3, # Determines the window length of smoothing that the savgol-filter uses for smoothing
        'smooth_polyorder': 2, # Determines the order of the polynomial used in the smoothing algorithm
        'smooth_algorithm': 'savgol', # At the present, only Savitzky-Golay filter is implemented. Add Gaussian and Boxcar later.
        'smooth_save_default': False, # Toggles whether or not to run a separate smoothing using default values on / off
        'smooth_store_data': False, # Toggles storing data to data['xanes_data'] on / off
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    # Initialise new DataFrame with correct x-values
    df_smooth = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    # Do the same if smoothing with default values is toggled on
    if options['smooth_save_default']:
        df_smooth_default = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    if options['log']:
        aux.write_log(message='Starting smoothing procedure.', options=options)


    # Run in interactive mode if enabled
    if options['interactive']:
        data['xanes_data_backup'] = data['xanes_data'] # Backup the data
        options['interactive'] = False # Turn interactive mode off so that it is not called again within the interactive loop
        options['show_plots'] = True # Force plotting on as interactive mode is useless without it
        smoothing_interactive(data=data, options=options) # Call interactive version of the function
        return


    # FIXME Add other types of filters
    for i, filename in enumerate(data['path']):

        if options['smooth_algorithm'] == 'savgol':
            if options['log']:
                aux.write_log(message=f'Smoothing {filename} with algorithm: {options["smooth_algorithm"]} ({i+1}/{len(data["path"])})', options=options)

            # Apply savgol filter and add to DataFrame
            df_smooth.insert(1, filename, savgol_filter(data['xanes_data'][filename], options['smooth_window_length'], options['smooth_polyorder']))
        
        if options['smooth_save_default']:
            if options['smooth_algorithm'] == 'savgol':
                if options['log']:
                    aux.write_log(message=f'Smoothing {filename} using default parameters with algorithm: {options["smooth_algorithm"]} ({i+1}/{len(data["path"])})', options=options)
                df_smooth_default.insert(1, filename, savgol_filter(data['xanes_data'][filename], default_options['smooth_window_length'], default_options['smooth_polyorder']))
        

        # Make plots ...
        if options['save_plots'] or options['show_plots']:

            

            edge_pos = estimate_edge_position(data=data, options=options)
            step_length = data['xanes_data']['ZapEnergy'].iloc[1] - data['xanes_data']['ZapEnergy'].iloc[0]


            # ... if default smoothing is enabled. Only plotting +- 10 step sizes from the edge position
            if options['smooth_save_default']:
                fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,5))
                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-10*step_length) & (data['xanes_data']['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='black', ax=ax1, kind='scatter')
                df_smooth.loc[(df_smooth['ZapEnergy'] > edge_pos-10*step_length) & (df_smooth['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
                ax1.set_title(f'{os.path.basename(filename)} - Smooth', size=20)

                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-10*step_length) & (data['xanes_data']['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='black', ax=ax2,  kind='scatter')
                df_smooth_default.loc[(df_smooth_default['ZapEnergy'] > edge_pos-10*step_length) & (df_smooth_default['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
                ax2.set_title(f'{os.path.basename(filename)} - Smooth (default values)', size=20)

            # ... if only smoothing with user defined variables is enabled. Only plotting +- 10 step sizes from the edge position
            elif not options['smooth_save_default']:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))
                data['xanes_data'].plot(x='ZapEnergy', y=filename, ax=ax1, kind='scatter', c='black')
                df_smooth.plot(x='ZapEnergy', y=filename, ax=ax1, c='red')

                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-10*step_length) & (data['xanes_data']['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='black', ax=ax2,  kind='scatter')
                df_smooth.loc[(df_smooth['ZapEnergy'] > edge_pos-10*step_length) & (df_smooth['ZapEnergy'] < edge_pos+10*step_length)].plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
                
                ax1.set_title(f'{os.path.basename(filename)} - Smooth', size=20)
                ax2.set_title(f'{os.path.basename(filename)} - Smooth Edge Region', size=20)

            # Save plots
            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_smooth.png'
                plt.savefig(dst, transparent=False)
            
            # Close plots
            if not options['show_plots']:
                plt.close()
    
    if not options['smooth_save_default']:
        df_smooth_default = None

    if options['smooth_store_data']:
        data['xanes_data'] = df_smooth
        options['smooth_store_data'] = False
    
    return df_smooth, df_smooth_default



def smoothing_interactive(data: dict, options: dict) -> None:
    ''' Defines the widgets to use with the ipywidgets interactive mode and calls the update function found in btp.ipywidgets. '''

    w = widgets.interactive(
        btp.ipywidgets_update, func=widgets.fixed(smoothing), data=widgets.fixed(data), options=widgets.fixed(options), 
        smooth_window_length=widgets.IntSlider(value=options['smooth_window_length'], min=3, max=21, step=2),
        smooth_polyorder=widgets.IntSlider(value=options['smooth_polyorder'],  min=1, max=5, step=1),
        smooth_store_data=widgets.Checkbox(value=options['smooth_store_data'])
    )
    
    options['widget'] = w

    display(w)

def backup(data):

    data['xanes_data_backup'] = data['xanes_data'].copy()


def restore_from_backup(data):
    ''' Restores DataFrame from data['xanes_data_backup'] to data['xanes_data']. This can be useful e.g. when smoothing and you want to re-do the smoothing with different parameters. 
    
    If there is no DataFrame stored in data['xanes_data_backup'], this function does nothing. '''

    if 'xanes_data_bakcup' in data.keys():
        data['xanes_data'] = data['xanes_data_backup'].copy()


def find_nearest(array, value):
    ''' Finds the value closest to value in array'''

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def estimate_edge_position(data: dict, options={}, index=0):
    ''' Gets an estimation of the edge position. This is very similar to determine_edge_position, but provides instead a quick and dirty way where the actual data point closest to the maximum of the differentiated data
    is located. '''

    required_options = ['log','logfile', 'periods']
    default_options = {
        'log': False, # Toggles logging on/off
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_edge_position_estimation.log', # Sets path to log-file
        'periods': 2, #Periods needs to be an even number for the shifting of values to work properly
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    #making new dataframe to keep the differentiated data
    df_diff = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])
    df_diff[data['path'][index]]=data['xanes_data_original'][data['path'][index]].diff(periods=options['periods'])

    #shifting column values up so that average differential fits right between the points used in the calculation
    df_diff[data['path'][index]]=df_diff[data['path'][index]].shift(-int(options['periods']/2))
    

    if 'pre_edge_masks' in options.keys():
        for mask in options['pre_edge_masks']:
            df_diff[data['path'][index]].loc[(df_diff['ZapEnergy'] > mask[0]) & (df_diff['ZapEnergy'] < mask[1])] = 0

    if 'post_edge_masks' in options.keys():
        for mask in options['post_edge_masks']:
            df_diff[data['path'][index]].loc[(df_diff['ZapEnergy'] > mask[0]) & (df_diff['ZapEnergy'] < mask[1])] = 0

    if 'edge_masks' in options.keys():
        for mask in options['edge_masks']:
            df_diff[data['path'][index]].loc[(df_diff['ZapEnergy'] > mask[0]) & (df_diff['ZapEnergy'] < mask[1])] = 0

    df_diff_max = df_diff[data['path'][index]].dropna().max()

    estimated_edge_pos = df_diff.loc[df_diff[data['path'][index]] == df_diff_max,'ZapEnergy'].values[0]

    if options['log']:
        aux.write_log(message=f'Estimated edge position is: {estimated_edge_pos} keV', options=options)

    return estimated_edge_pos

def determine_edge_position(data: dict, options={}):
    ''' Determines the edge position by 1) first differential maximum and/or 2) second differential zero-point. Calculates differential and/or double differential by diff.periods and double_diff.periods respectively.
    The differentiated and/or doubly differentiated data is fitted to a polynomial of diff.polyorder and/or double_diff.polyorder around the estimated edge position. The estimated edge position is set to be the x-value of the data 
    point at maximum of the differentiated data. The region to be fitted to the polynomial is determined by fit_region, which defaults to 5 times the distance between two data points, giving five data points to fit to.
    
    Allows plotting and saving of three plots to assess the quality of the fit, and also allows logging.
    
    Requires that XANES-data is already loaded in data['xanes_data']. This allows the user to choose when to determine the edge position - whether before or after normalisation, flattening etc.'''
    
    required_options = ['save_values', 'log', 'logfile', 'show_plots', 'save_plots', 'save_folder',  'diff', 'diff.polyorder', 'diff.periods', 'double_diff', 'double_diff.polyorder', 'double_diff.periods', 'points_around_edge', 'save_diff_data']
    default_options = {
        'save_values': True, # Whether the edge positions should be stored in a dictionary within the main data dictionary. 
        'log': False, # Toggles logging on/off
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_determine_edge_position.log', # Sets the path to the logfile. Ignored if log == False
        'show_plots': False, #  Toggles on/off whether plots should be shown. For sequential data, saving the plots and inspecting them there is probably better.
        'save_plots': False, # Toggles on/off whether plots should be saved. 
        'save_folder': './', # Sets the path to where the plots should be saved. Creates folder if doesn't exist. Ignored if save_plots == False
        'edge_masks': [],
        'diff': True, # Toggles calculation of the edge position based on differential data
        'diff.polyorder': 2, # Sets the order of the polynomial to fit edge region of the differential to
        'diff.periods': 2, # Sets the number of data points between which the first order difference should be calculated. Needs to be even for subsequent shifting of data to function.
        'double_diff': False, # Toggles calculation of the edge position based on double differential data
        'double_diff.polyorder': 1, # Sets the order of the polynomial to fit edge region of the double differential to
        'double_diff.periods': 2, # Sets the number of data points between which the second order difference should be calculated. Needs to be even for subsequent shifting of data to function.
        'points_around_edge': 1, # The length of the region to find points to fit to a function
        'save_diff_data': False
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    # Check if periods are even
    if options['diff'] and options['diff.periods'] % 2 != 0:
        if options['log']:
            aux.write_log(message='Periods for differentiation is not even. Ending run.', options=options)
        raise Exception("NB! Periods needs to be an even number for the shifting of values to work properly")
    if options['double_diff'] and options['double_diff.periods'] % 2 != 0:
        aux.write_log(message='Periods for double differentiation is not even. Ending run.', options=options)
        raise Exception("NB! Periods needs to be an even number for the shifting of values to work properly")


    if options['interactive']:
        data['xanes_data_backup'] = data['xanes_data']
        options['interactive'] = False
        options['interactive_session_active'] = True
        options['show_plots'] = True
        determine_edge_position_interactive(data=data, options=options)
        return


   
    
    # Prepare dataframes for differential data
    if options['diff']:
        df_diff = pd.DataFrame(data['xanes_data']['ZapEnergy'])
    if options['double_diff']:
        df_double_diff = pd.DataFrame(data['xanes_data']['ZapEnergy'])
    if options['save_values']:
        data['e0_diff'] = {}
        data['e0_double_diff'] = {}


    if options['log']:
        aux.write_log(message='Starting edge position determination', options=options)


    # Get rough estimate of edge position
    for i, filename in enumerate(data['path']):

        estimated_edge_pos = estimate_edge_position(data, options=options, index=i)

        
        fit_region = (options['points_around_edge']+1)*(data['xanes_data']['ZapEnergy'].iloc[1] - data['xanes_data']['ZapEnergy'].iloc[0])

        if fit_region < 0:
            fit_region = (options['points_around_edge']+1)*(data['xanes_data']['ZapEnergy'].iloc[10] - data['xanes_data']['ZapEnergy'].iloc[9])


        #========================== Fitting the first order derivative ==========

        if options['diff']:
            df_diff[filename] = data['xanes_data'][filename].diff(periods=options['diff.periods']) 
            df_diff[filename]=df_diff[filename].shift(-int(options['diff.periods']/2)) # Shifts the data back so that the difference between the points is located in the middle of the two points the caluclated difference is between

            # Picks out the points to be fitted
            df_diff_edge = df_diff.loc[(df_diff["ZapEnergy"] <= estimated_edge_pos+fit_region) & ((df_diff["ZapEnergy"] >= estimated_edge_pos-fit_region))]
    
            
            # Fitting a function to the chosen interval
            params = np.polyfit(df_diff_edge["ZapEnergy"], df_diff_edge[filename], options['diff.polyorder'])
            diff_function = np.poly1d(params)

            x_diff=np.linspace(df_diff_edge["ZapEnergy"].iloc[0],df_diff_edge["ZapEnergy"].iloc[-1],num=10000)
            y_diff=diff_function(x_diff)

            df_diff_fit_function = pd.DataFrame(x_diff)
            df_diff_fit_function['y_diff'] = y_diff
            df_diff_fit_function.columns = ['x_diff', 'y_diff']
            
            # Picks out the x-value where the y-value is at a maximum
            edge_pos_diff=x_diff[np.where(y_diff == np.amax(y_diff))][0]
            
            if options['log']:
                aux.write_log(message=f"... Edge position of {os.path.basename(filename)} determined by the differential maximum is: {str(round(edge_pos_diff,5))} keV", options=options)
            
            if options['save_values']:
                data['e0_diff'][filename] = edge_pos_diff

         #========================== Fitting the second order derivative ==========
        if options['double_diff']:
            df_double_diff[filename] = data['xanes_data'][filename].diff(periods=options['double_diff.periods']).diff(periods=options['double_diff.periods'])
            df_double_diff[filename]=df_double_diff[filename].shift(-int(options['double_diff.periods']))
            
            # Pick out region of interest
            df_double_diff_edge = df_double_diff.loc[(df_double_diff["ZapEnergy"] < estimated_edge_pos+fit_region) & ((df_double_diff["ZapEnergy"] > estimated_edge_pos-fit_region))]

            # Fitting a function to the chosen interval
            params = np.polyfit(df_double_diff_edge["ZapEnergy"], df_double_diff_edge[filename], options['double_diff.polyorder'])
            double_diff_function = np.poly1d(params)

            x_double_diff=np.linspace(df_double_diff_edge["ZapEnergy"].iloc[0], df_double_diff_edge["ZapEnergy"].iloc[-1],num=10000)
            y_double_diff=double_diff_function(x_double_diff)

            df_double_diff_fit_function = pd.DataFrame(x_double_diff)
            df_double_diff_fit_function['y_diff'] = y_double_diff
            df_double_diff_fit_function.columns = ['x_diff', 'y_diff']


            # Picks out the x-value where the y-value is closest to 0
            edge_pos_double_diff=x_double_diff[np.where(y_double_diff == find_nearest(y_double_diff,0))][0]
        
            if options['log']:
                aux.write_log(message=f"... Edge position of {os.path.basename(filename)} determined by the double differential zero-point is {str(round(edge_pos_double_diff,5))} keV", options=options)

                if options['diff']:
                    aux.write_log(message=f"... Difference between edge position estimated from differential maximum and double differential zero-point is {(edge_pos_diff-edge_pos_double_diff)*1000} eV.", options=options)

            if options['save_values']:
                data['e0_double_diff'][filename] = edge_pos_double_diff


        # Make and show / save plots ...
        if options['save_plots'] or options['show_plots']:


            # ... if both are enabled
            if options['diff'] and options['double_diff']:

                _, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(ncols=3, nrows=2, figsize=(20,20))
                data['xanes_data'].plot(x='ZapEnergy', y=filename, ax=ax1, c='black')
                ax1.axvline(x=edge_pos_diff, ls='--', c='green')
                
                df_diff.plot(x='ZapEnergy', y=filename, ax=ax2, kind='scatter')
                df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax2)
                ax2.set_xlim([edge_pos_diff-fit_region*1.5, edge_pos_diff+fit_region*1.5])
                ax2.axvline(x=estimated_edge_pos-fit_region, ls='--', c='black')
                ax2.axvline(x=edge_pos_diff, ls='--', c='green')
                ax2.axvline(x=estimated_edge_pos+fit_region, ls='--', c='black')
                ax2.set_title('Fit region of differentiated data')

                df_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax3, kind='scatter')
                df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax3)
                ax3.axvline(x=edge_pos_diff, ls='--', c='green')
                ax3.axvline(x=estimated_edge_pos, ls='--', c='red')
                ax3.set_title('Fit of differentiated data')


                data['xanes_data'].plot(x='ZapEnergy', y=filename, ax=ax4, c='black')
                ax4.axvline(x=edge_pos_double_diff, ls='--', c='green')

                df_double_diff.plot(x='ZapEnergy', y=filename, ax=ax5, kind='scatter')
                df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax5)
                ax5.set_xlim([edge_pos_double_diff-0.0015, edge_pos_double_diff+0.0015])
                ax5.axvline(x=estimated_edge_pos-fit_region, ls='--', c='black')
                ax5.axvline(x=edge_pos_double_diff, ls='--', c='green')
                ax5.axvline(x=estimated_edge_pos+fit_region, ls='--', c='black')

                df_double_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax6, kind='scatter')
                df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax6)
                ax6.axvline(x=edge_pos_double_diff, ls='--', c='green')
                ax6.axvline(x=estimated_edge_pos, ls='--', c='red')
                

            # ... if only first order differentials is enabled
            elif options['diff']:
                _, (ax1, ax2, ax3) = plt.subplots(ncols=3,nrows=1, figsize=(20, 10))
                
                data['xanes_data'].plot(x='ZapEnergy', y=filename, ax=ax1, c='black')
                ax1.axvline(x=edge_pos_diff, ls='--', c='green')

                df_diff.plot(x='ZapEnergy', y=filename, ax=ax2, kind='scatter')
                df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax2)
                ax2.set_xlim([edge_pos_diff-fit_region*1.5, edge_pos_diff+fit_region*1.5])
                ax2.axvline(x=edge_pos_diff-fit_region, ls='--', c='black')
                ax2.axvline(x=edge_pos_diff, ls='--', c='green')
                ax2.axvline(x=edge_pos_diff+fit_region, ls='--', c='black')

                df_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax3)
                df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax3)
                ax3.axvline(x=edge_pos_diff, ls='--', c='green')
                ax3.axvline(x=estimated_edge_pos, ls='--', c='red')

            # ... if only second order differentials is enabled
            elif options['double_diff']:
                _, (ax1, ax2, ax3) = plt.subplots(ncols=3,nrows=1, figsize=(20, 10))
                
                data['xanes_data'].plot(x='ZapEnergy', y=filename, ax=ax1, c='black')
                ax1.axvline(x=edge_pos_double_diff, ls='--', c='green')

                df_double_diff.plot(x='ZapEnergy', y=filename, ax=ax2, kind='scatter')
                df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax2)
                ax2.set_xlim([edge_pos_double_diff-fit_region*1.5, edge_pos_double_diff+fit_region*1.5])
                ax2.axvline(x=edge_pos_double_diff-fit_region, ls='--', c='black')
                ax2.axvline(x=edge_pos_double_diff, ls='--', c='green')
                ax2.axvline(x=edge_pos_double_diff+fit_region, ls='--', c='black')

                df_double_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax3)
                df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax3)
                ax3.axvline(x=edge_pos_double_diff, ls='--', c='green')
                ax3.axvline(x=estimated_edge_pos, ls='--', c='red')


            # Save plots if toggled
            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_edge_position.png'

                plt.savefig(dst, transparent=False)


            # Close plots if show_plots not toggled
            if not options['show_plots']:
                plt.close()


    if not options['diff']:
        edge_pos_diff = None
    if not options['double_diff']:
        edge_pos_double_diff = None


    if options['save_diff_data']:
        data['diff_data'] = df_diff if options['diff'] else None
        data['double_diff_data'] = df_double_diff if options['double_diff'] else None

    return edge_pos_diff, edge_pos_double_diff



def determine_edge_position_interactive(data: dict, options: dict) -> None:
    ''' Defines the widgets to use with the ipywidgets interactive mode and calls the update function found in btp.ipywidgets. '''

    w = widgets.interactive(
        btp.ipywidgets_update, func=widgets.fixed(determine_edge_position), data=widgets.fixed(data), options=widgets.fixed(options), 
        points_around_edge=widgets.IntSlider(value=options['points_around_edge'], min=1, max=20, step=1),   
    )
    
    options['widget'] = w

    display(w)

def determine_edge_shift(data: dict, options: dict, edge_pos: float) -> None:
    
    if 'edge' not in data.keys():
        data['edge'] = find_element(data)

    
    reference_energy = xas.edges.K['keV'].loc[xas.edges.K['Atom'] == data['edge']].values[0]

    edge_shift = reference_energy - edge_pos

    if options['log']:
        aux.write_log(message=f'Edge shift vs. reference value for {data["edge"]} is {edge_shift*1000} eV', options=options)

    return edge_shift

def normalise(data: dict, options={}):
    ''' Normalises the data so that the difference between the fitted pre- and post-edge functions is 1 at the edge position. 
    
    Requires that edge positions have already been determined with determine_edge_position() and stored in data['e0_diff']. '''


    required_options = ['log', 'logfile', 'normalisation_store_data']
    default_options = {
        'log': False, # Toggles logging on/off
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_normalisation.log', # Sets path to log-file
        'show_plots': False, #  Toggles on/off whether plots should be shown. For sequential data, saving the plots and inspecting them there is probably better.
        'save_plots': False, # Toggles on/off whether plots should be saved. 
        'save_folder': './', # Sets the path to where the plots should be saved. Creates folder if doesn't exist. Ignored if save_plots == False
        'normalisation_store_data': False, # Toggles storing of the flattened data in data['xanes_data'] on/off
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    normalised_df = pd.DataFrame(data['xanes_data']['ZapEnergy'])
    data['normalisation_constants'] = {}

    if options['normalisation_store_data']:
        pre_edge_fit_data_norm = pd.DataFrame(data['pre_edge_fit_data']['ZapEnergy'])
        post_edge_fit_data_norm = pd.DataFrame(data['post_edge_fit_data']['ZapEnergy'])

    #Finding the normalisation constant µ_0(E_0), by subtracting the value of the pre-edge-line from the value of the post-edge line at e0
    for filename in data['path']:
        e0_ind = data['post_edge_fit_data'].loc[data['post_edge_fit_data']['ZapEnergy'] == find_nearest(data['post_edge_fit_data']['ZapEnergy'], data['e0_diff'][filename])].index.values[0]
       
        #norm = data['post_edge_fit_data'][filename].iloc[find_nearest(data['post_edge_fit_data'][filename], data['e0'][filename])]
        normalisation_constant = data['post_edge_fit_data'][filename].iloc[e0_ind] - data['pre_edge_fit_data'][filename].iloc[e0_ind]
        print(normalisation_constant)
        normalised_df.insert(1, filename, data['xanes_data'][filename] / normalisation_constant)

    
        if options['show_plots'] or options['save_plots']:

            fig, ax = plt.subplots(figsize=(10,5))

            normalised_df.plot(x='ZapEnergy', y=filename, ax=ax, color='red', label='Normalised data')
            ax.set_title(f'{os.path.basename(filename)} - After normalisation', size=20)
            ax.set_ylabel('Normalised x$\mu$(E)', size=20)
            ax.set_xlabel('Energy (keV)', size=20)
            ax.axhline(y=1, ls='--', c='black')


            # Save plots if toggled
            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_normalisation.png'

                plt.savefig(dst, transparent=False)


            # Close plots if show_plots not toggled
            if not options['show_plots']:
                plt.close()


        if options['normalisation_store_data']:
            pre_edge_fit_data_norm.insert(1, filename, data['pre_edge_fit_data'][filename] / normalisation_constant)
            post_edge_fit_data_norm.insert(1, filename, data['post_edge_fit_data'][filename] / normalisation_constant)





    if options['normalisation_store_data']:
        data['xanes_data'] = normalised_df
        # Normalise the pre-edge and post-edge fit function data
        data['pre_edge_fit_data_norm'] = pre_edge_fit_data_norm
        data['post_edge_fit_data_norm'] = post_edge_fit_data_norm

        data['normalisation_constants'][filename] = normalisation_constant


    return normalised_df

    
def flatten(data:dict, options={}):
    ''' Flattens the post-edge region (from edge position and up). Only for visual purposes.
    
    Requires data['xanes_data'] that is normalised through normalise() and that normalised versions of the post_edge_fit_data is stored in data['post_edge_fit_data_norm'].
    Also assumes that the pre edge-fit data is already subtracted from the data'''


    required_options = ['log', 'logfile', 'flatten_store_data']
    default_options = {
        'log': False, # Toggles logging on/off
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_flattening.log', # Sets path to log-file
        'flatten_store_data': False, # Toggles storing of the flattened data in data['xanes_data'] on/off 
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    # Initialise DataFrame with x-values
    flattened_df = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    # Loop through all files
    for filename in data['path']:

        # Subtract 1 from the _normalised_ post edge fit function
        fit_function_diff = data['post_edge_fit_data_norm'][filename] - 1 - data['pre_edge_fit_data_norm'][filename]
        
        # Set all values from edge position and downwards to 0 so that only data above the edge position will be adjusted
        fit_function_diff.loc[flattened_df['ZapEnergy'] <= data['e0_diff'][filename]] = 0

        # Subtract the difference between 1 and the post edge fit function from the normalised data.
        flattened_df[filename] = data['xanes_data'][filename] - fit_function_diff


        if options['show_plots'] or options['save_plots']:

            fig, ax = plt.subplots(figsize=(10,5))

            flattened_df.plot(x='ZapEnergy', y=filename, ax=ax, color='red', label='Flattened data')
            ax.set_title(f'{os.path.basename(filename)} - After flattening', size=20)
            ax.set_ylabel('Normalised x$\mu$(E)', size=20)
            ax.set_xlabel('Energy (keV)', size=20)
            ax.axhline(y=1, ls='--', c='black')


            # Save plots if toggled
            if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_flattened.png'

                plt.savefig(dst, transparent=False)


            # Close plots if show_plots not toggled
            if not options['show_plots']:
                plt.close()


    # Saves the flattened DataFrame
    if options['flatten_store_data']:
        data['xanes_data'] = flattened_df
    

    return flattened_df, fit_function_diff




