from logging import raiseExceptions
from jinja2 import TemplateRuntimeError
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import nafuma.auxillary as aux
import nafuma.xanes as xas
import nafuma.xanes.io as io
from scipy.signal import savgol_filter
from datetime import datetime


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

    required_options = ['pre_edge_start', 'log', 'logfile', 'save_plots', 'save_folder']
    default_options = {
        'pre_edge_start': None,
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_pre_edge_fit.log',
        'save_plots': False,
        'save_folder': './'
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if options['log']:
        aux.write_log(message='Starting pre edge fit', options=options)



    # FIXME Implement with finding accurate edge position
    # FIXME Allow specification of start of pre-edge area
    # Find the cutoff point at which the edge starts - everything to the LEFT of this point will be used in the pre edge function fit
    if not options['pre_edge_start']:
        pre_edge_limit_offset = 0.03

        data['edge'] = find_element(data)

        edge_position = estimate_edge_position(data, options, index=0)
        pre_edge_limit = edge_position - pre_edge_limit_offset

    # FIXME There should be an option to specify the interval in which to fit the background - now it is taking everything to the left of edge_start parameter, but if there are some artifacts in this area, it should be possible to
    # limit the interval
    # Making a dataframe only containing the rows that are included in the background subtraction (points lower than where the edge start is defined)
    pre_edge_data = data['xanes_data_original'].loc[data['xanes_data_original']["ZapEnergy"] < pre_edge_limit]
        
    # Making a new dataframe, with only the ZapEnergies as the first column -> will be filled to include the background data
    pre_edge_fit_data = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    for i, filename in enumerate(data['path']):
        if options['log']:
            aux.write_log(message=f'Fitting background on {os.path.basename(filename)} ({i+1} / {len(data["path"])})', options=options)

        #Fitting linear function to the background
        params = np.polyfit(pre_edge_data["ZapEnergy"],pre_edge_data[filename],1)
        fit_function = np.poly1d(params)
        
        #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        background=fit_function(pre_edge_fit_data["ZapEnergy"])
            
        #adding a new column in df_background with the y-values of the background
        pre_edge_fit_data.insert(1,filename,background) 
        
        if options['save_plots']:
            if not os.path.isdir(options['save_folder']):
                os.makedirs(options['save_folder'])

            dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_pre_edge_fit.png'

            fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax1)
            pre_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
            ax1.axvline(x = max(pre_edge_data['ZapEnergy']), ls='--')
            ax1.set_title(f'{os.path.basename(filename)} - Full view', size=20)

            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax2)
            pre_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
            ax2.axvline(x = max(pre_edge_data['ZapEnergy']), ls='--')
            ax2.set_xlim([min(pre_edge_data['ZapEnergy']), max(pre_edge_data['ZapEnergy'])])
            ax2.set_ylim([min(pre_edge_data[filename]), max(pre_edge_data[filename])])
            ax2.set_title(f'{os.path.basename(filename)} - Fit region', size=20)


            plt.savefig(dst, transparent=False)
            plt.close()


    if options['log']:
        aux.write_log(message=f'Pre edge fitting done.', options=options)   

    return pre_edge_fit_data



def pre_edge_subtraction(data: dict, options={}):

    required_options = ['log', 'logfile', 'save_plots', 'save_folder']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_pre_edge_subtraction.log',
        'save_plots': False,
        'save_folder': './'
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    if options['log']:
        aux.write_log(message='Starting pre edge subtraction', options=options)

    xanes_data_bkgd_subtracted = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    for i, filename in enumerate(data['path']):
        if options['log']:
            aux.write_log(message=f'Subtracting background on {filename} ({i} / {len(data["path"])}', options=options)

        xanes_data_bkgd_subtracted.insert(1, filename, data['xanes_data_original'][filename] - data['pre_edge_fit_data'][filename])

        if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_pre_edge_subtraction.png'

                fig, ax = plt.subplots(figsize=(10,5))
                data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax)
                xanes_data_bkgd_subtracted.plot(x='ZapEnergy', y=filename, color='red', ax=ax)
                ax.set_title(f'{os.path.basename(filename)} - After subtraction', size=20)

                plt.savefig(dst)
                plt.close()

    return xanes_data_bkgd_subtracted





def post_edge_fit(data: dict, options={}):
    #FIXME should be called "fitting post edge" (normalization is not done here, need edge shift position)
    required_options = ['log', 'logfile', 'post_edge_interval']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_post_edge_fit.log',
        'post_edge_interval': [None, None],
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if not options['post_edge_interval'][0]:
        post_edge_limit_offset = 0.03

        data['edge'] = find_element(data)

        edge_position = estimate_edge_position(data, options, index=0)
        options['post_edge_interval'][0] = edge_position + post_edge_limit_offset


    if not options['post_edge_interval'][1]:
        options['post_edge_interval'][1] = data['xanes_data_original']['ZapEnergy'].max()


    post_edge_data = data['xanes_data_original'].loc[(data['xanes_data_original']["ZapEnergy"] > options['post_edge_interval'][0]) & (data['xanes_data_original']["ZapEnergy"] < options['post_edge_interval'][1])]
    post_edge_data.dropna(inplace=True) #Removing all indexes without any value, as some of the data sets misses the few last data points and fucks up the fit

    # Making a new dataframe, with only the ZapEnergies as the first column -> will be filled to include the background data
    post_edge_fit_data = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])
    
    for i, filename in enumerate(data['path']):
        if options['log']:
            aux.write_log(message=f'Fitting post edge on {os.path.basename(filename)} ({i+1} / {len(data["path"])})', options=options)

        #Fitting linear function to the background
        params = np.polyfit(post_edge_data["ZapEnergy"], post_edge_data[filename], 2)
        fit_function = np.poly1d(params)
        
        #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        background=fit_function(post_edge_fit_data["ZapEnergy"])
            
        #adding a new column in df_background with the y-values of the background
        post_edge_fit_data.insert(1,filename,background) 
        
        if options['save_plots']:
            if not os.path.isdir(options['save_folder']):
                os.makedirs(options['save_folder'])

            dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_post_edge_fit.png'

            fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax1)
            post_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
            ax1.axvline(x = max(post_edge_data['ZapEnergy']), ls='--')
            ax1.set_title(f'{os.path.basename(filename)} - Full view', size=20)

            data['xanes_data_original'].plot(x='ZapEnergy', y=filename, color='black', ax=ax2)
            post_edge_fit_data.plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
            ax2.axvline(x = max(post_edge_data['ZapEnergy']), ls='--')
            ax2.set_xlim([min(post_edge_data['ZapEnergy']), max(post_edge_data['ZapEnergy'])])
            ax2.set_ylim([min(post_edge_data[filename]), max(post_edge_data[filename])])
            ax2.set_title(f'{os.path.basename(filename)} - Fit region', size=20)


            plt.savefig(dst, transparent=False)
            plt.close()


    return post_edge_fit_data

def smoothing(data: dict, options={}):

    # FIXME Add logging
    # FIXME Add saving of files

    required_options = ['log', 'logfile', 'window_length','polyorder', 'save_default']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_smoothing.log',
        'save_plots': False,
        'save_folder': './',
        'window_length': 3,
        'polyorder': 2,
        'save_default': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    df_smooth = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    if options['save_default']:
        df_smooth_default = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    # FIXME Add other types of filters
    # FIXME Instead of assigning values directly to the data dictionary, these should be made into an own DataFrame that you can decide later what to do with - these variables should
    # then be returned
    for filename in data['path']:
        df_smooth.insert(1, filename, savgol_filter(data['xanes_data'][filename], options['window_length'], options['polyorder']))
        
        if options['save_default']:
            df_smooth_default.insert(1, filename, savgol_filter(data['xanes_data'][filename], default_options['window_length'], default_options['polyorder']))
        

        if options['save_plots']:
            if not os.path.isdir(options['save_folder']):
                os.makedirs(options['save_folder'])

            dst = os.path.join(options['save_folder'], os.path.basename(filename)) + '_smooth.png'

            edge_pos = estimate_edge_position(data=data, options=options)
            intensity_midpoint = df_smooth[filename].iloc[np.where(df_smooth['ZapEnergy'] == find_nearest(df_smooth['ZapEnergy'], edge_pos))].values[0]
    
            if options['save_default']:
                fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,5))
                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-0.0015) & (data['xanes_data']['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='black', ax=ax1, kind='scatter')
                df_smooth.loc[(df_smooth['ZapEnergy'] > edge_pos-0.0015) & (df_smooth['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='red', ax=ax1)
                ax1.set_title(f'{os.path.basename(filename)} - Smooth', size=20)

                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-0.0015) & (data['xanes_data']['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='black', ax=ax2,  kind='scatter')
                df_smooth_default.loc[(df_smooth_default['ZapEnergy'] > edge_pos-0.0015) & (df_smooth_default['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='red', ax=ax2)
                ax2.set_title(f'{os.path.basename(filename)} - Smooth (default values)', size=20)

            elif not options['save_default']:
                fig, ax = plt.subplots(figsize=(10,5))
                data['xanes_data'].loc[(data['xanes_data']['ZapEnergy'] > edge_pos-0.0015) & (data['xanes_data']['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='black', ax=ax,  kind='scatter')
                df_smooth.loc[(df_smooth['ZapEnergy'] > edge_pos-0.0015) & (df_smooth['ZapEnergy'] < edge_pos+0.0015)].plot(x='ZapEnergy', y=filename, color='red', ax=ax)
                ax.set_xlim([edge_pos-0.0015, edge_pos+0.0015])
                ax.set_ylim([intensity_midpoint*0.9, intensity_midpoint*1.1])
                
                ax.set_title(f'{os.path.basename(filename)} - Smooth', size=20)


            plt.savefig(dst, transparent=False)
            plt.close()
    
    if not options['save_default']:
        df_smooth_default = None
    
    return df_smooth, df_smooth_default



def find_nearest(array, value):
    #function to find the value closes to "value" in an "array"
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def estimate_edge_position(data: dict, options={}, index=0):
    #a dataset is differentiated to find a first estimate of the edge shift to use as starting point. 
    required_options = ['log','logfile', 'periods']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_edge_position_estimation.log',
        'periods': 2, #Periods needs to be an even number for the shifting of values to work properly
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    #making new dataframe to keep the differentiated data
    df_diff = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])
    df_diff[data['path'][index]]=data['xanes_data_original'][data['path'][index]].diff(periods=options['periods'])

    #shifting column values up so that average differential fits right between the points used in the calculation
    df_diff[data['path'][index]]=df_diff[data['path'][index]].shift(-int(options['periods']/2))
    df_diff_max = df_diff[data['path'][index]].dropna().max()
    estimated_edge_shift =df_diff.loc[df_diff[data['path'][index]] == df_diff_max,'ZapEnergy'].values[0]

    # FIXME Add logging option to see the result

    if options['log']:
        aux.write_log(message=f'Estimated edge shift for determination of pre-edge area is: {estimated_edge_shift} keV', options=options)

    return estimated_edge_shift

def determine_edge_position(data: dict, options={}):
    
    required_options = ['log', 'logfile', 'save_plots', 'save_folder', 'periods', 'diff', 'double_diff', 'fit_region']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_determine_edge_position.log',
        'save_plots': False,
        'save_folder': './',
        'periods': 2, #Periods needs to be an even number for the shifting of values to work properly,
        'diff': True,
        'double_diff': False,
        'fit_region': 0.0005
    
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if options['periods'] % 2 == 1:
        raise Exception("NB! Periods needs to be an even number for the shifting of values to work properly")

   
    #####

    if options['diff']:
        df_diff = pd.DataFrame(data['xanes_data']['ZapEnergy'])
    if options['double_diff']:
        df_double_diff = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    for i, filename in enumerate(data['path']):
        estimated_edge_pos = estimate_edge_position(data, options=options, index=i)

        
        #========================== fitting first differential ==========

        if options['diff']:
            df_diff[filename] = data['xanes_data'][filename].diff(periods=options['periods'])
            df_diff[filename]=df_diff[filename].shift(-int(options['periods']/2))

            df_diff_edge = df_diff.loc[(df_diff["ZapEnergy"] < estimated_edge_pos+options['fit_region']) & ((df_diff["ZapEnergy"] > estimated_edge_pos-options['fit_region']))]
    
            
            # Fitting a function to the chosen interval
            params = np.polyfit(df_diff_edge["ZapEnergy"], df_diff_edge[filename], 2)
            diff_function = np.poly1d(params)

            x_diff=np.linspace(df_diff_edge["ZapEnergy"].iloc[0],df_diff_edge["ZapEnergy"].iloc[-1],num=10000)
            y_diff=diff_function(x_diff)

            df_diff_fit_function = pd.DataFrame(x_diff)
            df_diff_fit_function['y_diff'] = y_diff
            df_diff_fit_function.columns = ['x_diff', 'y_diff']
            
            # Picks out the x-value where the y-value is at a maximum
            edge_pos_diff=x_diff[np.where(y_diff == np.amax(y_diff))][0]
            
            if options['log']:
                aux.write_log(message=f"Edge shift estimated by the differential maximum is: {str(round(edge_pos_diff,5))}", options=options)


        if options['double_diff']:
            df_double_diff[filename] = data['xanes_data'][filename].diff(periods=options['periods']).diff(periods=options['periods'])
            df_double_diff[filename]=df_double_diff[filename].shift(-int(options['periods']))
            
            # Pick out region of interest
            df_double_diff_edge = df_double_diff.loc[(df_double_diff["ZapEnergy"] < estimated_edge_pos+options['fit_region']) & ((df_double_diff["ZapEnergy"] > estimated_edge_pos-options['fit_region']))]

            # Fitting a function to the chosen interval
            params = np.polyfit(df_double_diff_edge["ZapEnergy"], df_double_diff_edge[filename], 2)
            double_diff_function = np.poly1d(params)

            x_double_diff=np.linspace(df_double_diff_edge["ZapEnergy"].iloc[0], df_double_diff_edge["ZapEnergy"].iloc[-1],num=10000)
            y_double_diff=double_diff_function(x_double_diff)

            df_double_diff_fit_function = pd.DataFrame(x_double_diff)
            df_double_diff_fit_function['y_diff'] = y_double_diff
            df_double_diff_fit_function.columns = ['x_diff', 'y_diff']


            # Picks out the x-value where the y-value is closest to 0
            edge_pos_double_diff=x_double_diff[np.where(y_double_diff == find_nearest(y_double_diff,0))][0]
        
            if options['log']:
                aux.write_log(message=f"Edge shift estimated by the double differential zero-point is {str(round(edge_pos_double_diff,5))}", options=options)

            if options['save_plots']:

                if options['diff'] and options['double_diff']:
                    
                    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(20,20))
                    df_diff.plot(x='ZapEnergy', y=filename, ax=ax1, kind='scatter')
                    df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax1)
                    ax1.set_xlim([edge_pos_diff-0.0015, edge_pos_diff+0.0015])
                    ax1.axvline(x=edge_pos_diff-options['fit_region'], ls='--', c='black')
                    ax1.axvline(x=edge_pos_diff, ls='--', c='green')
                    ax1.axvline(x=edge_pos_diff+options['fit_region'], ls='--', c='black')
                    ax1.set_title('Fit region of differentiated data')

                    df_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax2, kind='scatter')
                    df_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax2)
                    ax2.axvline(x=edge_pos_diff, ls='--', c='green')
                    ax2.axvline(x=estimated_edge_pos, ls='--', c='red')
                    ax2.set_title('Fit of differentiated data')


                    df_double_diff.plot(x='ZapEnergy', y=filename, ax=ax3, kind='scatter')
                    df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax3)
                    ax3.set_xlim([edge_pos_double_diff-0.0015, edge_pos_double_diff+0.0015])
                    ax3.axvline(x=edge_pos_double_diff-options['fit_region'], ls='--', c='black')
                    ax3.axvline(x=edge_pos_double_diff, ls='--', c='green')
                    ax3.axvline(x=edge_pos_double_diff+options['fit_region'], ls='--', c='black')

                    df_double_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax4, kind='scatter')
                    df_double_diff_fit_function.plot(x='x_diff', y='y_diff', ax=ax4)
                    ax4.axvline(x=edge_pos_double_diff, ls='--', c='green')
                    ax4.axvline(x=estimated_edge_pos, ls='--', c='red')

                    


                elif options['diff']:
                    fig, (ax1, ax2) = plt.subplots(ncols=2,nrows=1, figsize=(20, 10))
                    df_diff.plot(x='ZapEnergy', y=filename, ax=ax1, kind='scatter')
                    ax1.set_xlim([edge_pos_diff-0.5, edge_pos_diff+0.5])
                    ax1.axvline(x=edge_pos_diff-options['fit_region'], ls='--', c='black')
                    ax1.axvline(x=edge_pos_diff, ls='--', c='green')
                    ax1.axvline(x=edge_pos_diff+options['fit_region'], ls='--', c='black')

                    df_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax2)
                    ax2.axvline(x=edge_pos_diff, ls='--', c='green')
                    ax2.axvline(x=estimated_edge_pos, ls='--', c='red')

                
                elif options['double_diff']:
                    fig, (ax1, ax2) = plt.subplots(ncols=2,nrows=1, figsize=(20, 10))
                    df_double_diff.plot(x='ZapEnergy', y=filename, ax=ax1, kind='scatter')
                    ax1.set_xlim([edge_pos_double_diff-0.5, edge_pos_double_diff+0.5])
                    ax1.axvline(x=edge_pos_double_diff-options['fit_region'], ls='--', c='black')
                    ax1.axvline(x=edge_pos_double_diff, ls='--', c='green')
                    ax1.axvline(x=edge_pos_double_diff+options['fit_region'], ls='--', c='black')

                    df_double_diff_edge.plot(x='ZapEnergy', y=filename, ax=ax2)
                    ax2.axvline(x=edge_pos_double_diff, ls='--', c='green')
                    ax2.axvline(x=estimated_edge_pos, ls='--', c='red')


    if not options['diff']:
        edge_pos_diff = None
    if not options['double_diff']:
        edge_pos_double_diff = None

    return edge_pos_diff, edge_pos_double_diff

def normalise(data: dict, options={}):
    required_options = ['log', 'logfile', 'save_values']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_normalisation.log',
        'save_values': True
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    normalised_df = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    #Finding the normalisation constant µ_0(E_0), by subtracting the value of the pre-edge-line from the value of the post-edge line at e0
    for filename in data['path']:
        normalisation_constant = data['post_edge_fit_function'][filename].loc[data['post_edge_fit_function']['ZapEnergy'] == data['e0'][filename]] - data['pre_edge_fit_function'].loc[data['pre_edge_fit_function']['ZapEnergy'] == data['e0'][filename]]

        normalised_df.insert(1, filename, data['xanes_data'] / normalisation_constant)

    if options['save_values']:
        data['xanes_data'] = normalised_df


    return normalised_df

    
def flatten(data:dict, options={}):
    #only picking out zapenergy-values higher than edge position (edge pos and below remains untouched)

    required_options = ['log', 'logfile', 'save_values']
    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_flattening.log',
        'save_values': True
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    df_e0_and_above=df.loc[df['ZapEnergy'] > edge_shift_diff]

    flattened_df = pd.DataFrame(data['xanes_data']['ZapEnergy'])

    for filename in data['path']:
        above_e0 = data['xanes_data'][filename].loc(data['xanes_data']['ZapEnergy'] > data['e0'][filename])
        flattened_data = data['post_edge_fit_function'][filename] - 



    flattened_data = post_edge_fit_function(df_e0_and_above['ZapEnergy']) - pre_edge_fit_function(df_e0_and_above['ZapEnergy'])

    #make a new dataframe with flattened values


