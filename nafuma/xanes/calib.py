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
def find_element(data: dict) -> str:
    ''' Takes the data dictionary and determines based on the start value of the ZapEnergy-column which element the edge is from.'''

    element_energy_intervals = {
        'Mn': [5.9, 6.5],
        'Fe': [7.0, 7.2],
        'Co': [7.6, 7.8],
        'Ni': [8.0, 8.6]
    }

    if (element_energy_intervals['Mn'][0] < data['xanes_data_original']["ZapEnergy"].iloc[0]) & (data['xanes_data_original']["ZapEnergy"].iloc[0] < element_energy_intervals['Mn'][1]):
        edge = 'Mn'
    elif (element_energy_intervals['Fe'][0] < data['xanes_data_original']["ZapEnergy"].iloc[0]) & (data['xanes_data_original']["ZapEnergy"].iloc[0] < element_energy_intervals['Fe'][1]):
        edge = 'Fe'
    elif (element_energy_intervals['Co'][0] < data['xanes_data_original']["ZapEnergy"].iloc[0]) & (data['xanes_data_original']["ZapEnergy"].iloc[0] < element_energy_intervals['Co'][1]):
        edge = 'Co'   
    elif (element_energy_intervals['Ni'][0] < data['xanes_data_original']["ZapEnergy"].iloc[0]) & (data['xanes_data_original']["ZapEnergy"].iloc[0] < element_energy_intervals['Ni'][1]):
        edge = 'Ni'
        
        
    return(edge)



def pre_edge_fit(data: dict, options={}) -> pd.DataFrame:


    # FIXME Add log-file

    required_options = ['edge_start', 'log', 'logfile', 'save_plots', 'save_folder']
    default_options = {
        'edge_start': None,
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_pre_edge_fit.log',
        'save_plots': False,
        'save_folder': './'
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if options['log']:
        aux.write_log(message='Starting pre edge fit', options=options)



    # FIXME Implement with finding accurate edge position
    # Find the cutoff point at which the edge starts - everything to the LEFT of this point will be used in the pre edge function fit
    if not options['edge_start']:
        edge_starts = {
            'Mn': 6.42,
            'Fe': 7.09,
            'Co': 7.705,
            'Ni': 8.3        
        }

        data['edge'] = find_element(data)
        edge_start = edge_starts[data['edge']]

    # FIXME There should be an option to specify the interval in which to fit the background - now it is taking everything to the left of edge_start parameter, but if there are some artifacts in this area, it should be possible to
    # limit the interval
    # Making a dataframe only containing the rows that are included in the background subtraction (points lower than where the edge start is defined)
    pre_edge_data = data['xanes_data_original'].loc[data['xanes_data_original']["ZapEnergy"] < edge_start]
        
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
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S.log")}_pre_edge_subtraction.log',
        'save_plots': False,
        'save_folder': './'
    }

    if options['log']:
        aux.write_log(message='Starting pre edge subtraction', options=options)

    xanes_data_bkgd_subtracted = pd.DataFrame(data['xanes_data_original']["ZapEnergy"])

    for i, filename in enumerate(data['path']):
        if options['log']:
            aux.write_log(message=f'Subtracting background on {filename} ({i} / {len(data["path"])}', options=options)

        xanes_data_bkgd_subtracted.insert(1, filename, data['xanes_data'][filename] - data['pre_edge_fit_data'][filename])

        if options['save_plots']:
                if not os.path.isdir(options['save_folder']):
                    os.makedirs(options['save_folder'])

                dst = os.path.join(options['save_folder'], filename) + '_pre_edge_subtraction.png'

                fig, ax = plt.subplots(1,2,figsize=(10,5))
                data['xanes_data'].plot(x='ZapEnergy', y=filename, color='black', ax=ax)
                xanes_data_bkgd_subtracted.plot(x='ZapEnergy', y=filename, color='red', ax=ax)
                ax.set_title(f'{os.path.basename(filename)} - After subtraction', size=20)

                plt.savefig(dst)
                plt.close()

    return xanes_data_bkgd_subtracted


def post_edge_fit(path, options={}):
    #FIXME should be called "fitting post edge" (normalization is not done here, need edge shift position)
    required_options = ['print']
    default_options = {
        'print': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    
    df_bkgd_sub,filenames,edge = pre_edge_subtraction(path, options=options)
    #Defining the end of the pre-edge-region for Mn/Ni, thus start of the edge
    #FIXME Use rought edge shift estimate, add X eV as first guess, have an option to adjust this value with widget
    if edge == 'Mn':
        edge_stop = 6.565
    if edge == 'Ni':
        edge_stop = 8.361

    df_end= df_bkgd_sub.loc[df_bkgd_sub["ZapEnergy"] > edge_stop] # new dataframe only containing the post edge, where a regression line will be calculated in the for-loop below
    df_end.dropna(inplace=True) #Removing all indexes without any value, as some of the data sets misses the few last data points and fucks up the fit
    df_postedge = pd.DataFrame(df_bkgd_sub["ZapEnergy"]) #making a new dataframe 

    function_post_list=[]
    for files in filenames: 
        d = np.polyfit(df_end["ZapEnergy"],df_end[files],1)
        function_post = np.poly1d(d)
        y_post=function_post(df_bkgd_sub["ZapEnergy"]) 
        function_post_list.append(function_post)
        df_postedge.insert(1,files,y_post) #adding a new column with the y-values of the fitted post edge

    #Plotting the background subtracted signal with the post-edge regression line and the start point for the linear regression line
    if options['print'] == True:
        ax = df_bkgd_sub.plot(x = "ZapEnergy",y=filenames) #defining x and y
        plt.axvline(x = min(df_end["ZapEnergy"])) 
        fig = plt.figure(figsize=(15,15))
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        ax = df_bkgd_sub.plot(x = "ZapEnergy",y=filenames, legend=False) #defining x and y
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        plt.axvline(x = min(df_end["ZapEnergy"])) 

    return df_bkgd_sub, df_postedge, filenames, edge

def smoothing(path, options={}):
    required_options = ['print','window_length','polyorder']
    default_options = {
        'print': False,
        'window_length': 3,
        'polyorder': 2
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    df_bkgd_sub, df_postedge, filenames, edge = post_edge_fit(path,options=options)
    #================= SMOOTHING
    df_smooth = pd.DataFrame(df_bkgd_sub["ZapEnergy"])
    df_default = pd.DataFrame(df_bkgd_sub["ZapEnergy"])
    #df_smooth[filenames] = df_bkgd_sub.iloc[:,2].rolling(window=rolling_av).mean()
    #df_smooth[filenames] = df_smooth[filenames].shift(-int((rolling_av)/2))
    for filename in filenames:
        x_smooth=savgol_filter(df_bkgd_sub[filename], options['window_length'],options['polyorder'])
        df_smooth[filename] = x_smooth
        x_default=savgol_filter(df_bkgd_sub[filename],default_options['window_length'],default_options['polyorder'])
        df_default[filename] = x_default
    
    
        
    #printing the smoothed curves vs data
    if options['print'] == True:

        ## ================================================
        #df_diff = pd.DataFrame(df_smooth["ZapEnergy"])
        #df_diff_estimated_max = df_diff[filenames].dropna().max()

    
        #estimated_edge_shift=df_diff.loc[df_diff[filenames] == df_diff_max,'ZapEnergy'].values[0]
        # ==========================================
        

        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))
        x_range_zoom=[6.54,6.55] #make into widget
        y_range_zoom=[20000,80000] #make into widget

        df_bkgd_sub.plot.scatter(x = "ZapEnergy",y=filenames, ax=ax1, color="Red") 
        df_smooth.plot(x = "ZapEnergy",y=filenames, ax=ax1, color="Blue")
        ax1.set_xlim(x_range_zoom)
        ax1.set_ylim(y_range_zoom)
        ax1.set_title("Smoothed curve (blue) vs data (red) used for further analysis")

        df_bkgd_sub.plot.scatter(x = "ZapEnergy",y=filenames, ax=ax2, color="Red") 
        df_default.plot(x = "ZapEnergy",y=filenames, ax=ax2, color="Green") 
        ax2.set_xlim(x_range_zoom)
        ax2.set_ylim(y_range_zoom)
        ax2.set_title("Smoothed curve (green) vs data (red) using default window_length and polyorder")
    
    return df_smooth, filenames

def find_pos_maxdiff(df, filenames,options={}):
    #a dataset is differentiated to find a first estimate of the edge shift to use as starting point. 
    required_options = ['print','periods']
    default_options = {
        'print': False,
        'periods': 2, #Periods needs to be an even number for the shifting of values to work properly
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    #making new dataframe to keep the differentiated data
    df_diff = pd.DataFrame(df["ZapEnergy"])
    df_diff[filenames]=df[filenames].diff(periods=options['periods'])

    #shifting column values up so that average differential fits right between the points used in the calculation
    df_diff[filenames]=df_diff[filenames].shift(-int(options['periods']/2))
    df_diff_max = df_diff[filenames].dropna().max()
    estimated_edge_shift =df_diff.loc[df_diff[filenames] == df_diff_max,'ZapEnergy'].values[0]

    return estimated_edge_shift, df_diff, df_diff_max

def find_nearest(array, value):
    #function to find the value closes to "value" in an "array"
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def finding_e0(path, options={}):
    required_options = ['print','periods']
    default_options = {
        'print': False,
        'periods': 2, #Periods needs to be an even number for the shifting of values to work properly
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    df_smooth, filenames = smoothing(path, options=options) #This way the smoothing is printed as long as the "finding e0" is printed.

    if options['periods'] % 2 == 1:
        print("NB!!!!!!!!!!!!!!!!! Periods needs to be an even number for the shifting of values to work properly")
    ###df_diff = pd.DataFrame(df_smooth["ZapEnergy"]) #
    if len(filenames) == 1:
        filenames=filenames[0]
    else:
       print("MORE THAN ONE FILE --> generalize")
    
    #####
    estimated_edge_shift, df_diff, df_diff_max = find_pos_maxdiff(df_smooth, filenames,options=options)
    print(estimated_edge_shift)
    ####
    ###df_diff[filenames]=df_smooth[filenames].diff(periods=options['periods']) #
    df_doublediff=pd.DataFrame(df_smooth["ZapEnergy"])
    df_doublediff[filenames]=df_diff[filenames].diff(periods=options['periods'])
    
    if options['print'] == True:
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))

        df_diff.plot(x = "ZapEnergy",y=filenames, ax=ax1) #defining x and y
        df_doublediff.plot(x = "ZapEnergy",y=filenames,ax=ax2) #defining x and y

    #shifting column values up so that average differential fits right between the points used in the calculation
    #df_diff[filenames]=df_diff[filenames].shift(-int(options['periods']/2)) #
    df_doublediff[filenames]=df_doublediff[filenames].shift(-int(options['periods']))
    
    #finding maximum value to maneuver to the correct part of the data set
    #df_diff_max = df_diff[filenames].dropna().max()

    
    estimated_edge_shift=df_diff.loc[df_diff[filenames] == df_diff_max,'ZapEnergy'].values[0]

    fit_region = 0.0004
    df_diff_edge=df_diff.loc[(df_diff["ZapEnergy"] < estimated_edge_shift+fit_region)]# and (df_diff["ZapEnergy"] > estimated_edge_shift-0.05)]
    df_diff_edge=df_diff_edge.loc[(df_diff["ZapEnergy"] > estimated_edge_shift-fit_region)]
    
    
    
    
    df_doublediff_edge=df_doublediff.loc[(df_doublediff["ZapEnergy"] < estimated_edge_shift+fit_region)]# and (df_diff["ZapEnergy"] > estimated_edge_shift-0.05)]
    df_doublediff_edge=df_doublediff_edge.loc[(df_doublediff["ZapEnergy"] > estimated_edge_shift-fit_region)]
    #df_diff_edge=df_diff.loc[(df_diff["ZapEnergy"] > estimated_edge_shift-0.15) and (df_diff["ZapEnergy"] < estimated_edge_shift+0.15)]

    #df_diff_edge=df_diff.loc[df_diff["ZapEnergy"] > estimated_edge_shift-0.15]
    #print(df_diff_edge)
    if options['print'] == True:
        fig, (ax3,ax4) = plt.subplots(1,2,figsize=(15,5))

        df_diff_edge.plot(x = "ZapEnergy",y=filenames,ax=ax3) #defining x and y
        ax3.set_title("Zoomed into edge region (derivative))")
        ax3.axvline(x = estimated_edge_shift)

        df_doublediff_edge.plot(x = "ZapEnergy",y=filenames,ax=ax4,kind="scatter") #defining x and y
        ax4.set_title("Zoomed into edge region (double derivative)")
        ax4.axvline(x = estimated_edge_shift)
        ax4.axhline(0)



        #ax1.set_xlim([estimated_edge_shift-fit_region,estimated_edge_shift+fit_region])
        #ax1.set_title("not sure what this is tbh")
  
        #ax2.set_xlim([estimated_edge_shift-fit_region,estimated_edge_shift+fit_region])
        #ax2.set_title("not sure what this is either tbh")

    #==============
    #df_smooth=df_smooth2
    #=================




    #========================== fitting first differential ==========
    df_diff = df_diff[df_diff[filenames].notna()]

    #fitting a function to the chosen interval
    d = np.polyfit(df_diff_edge["ZapEnergy"],df_diff_edge[filenames],2)
    function_diff = np.poly1d(d)

    x_diff=np.linspace(df_diff_edge["ZapEnergy"].iloc[0],df_diff_edge["ZapEnergy"].iloc[-1],num=1000)
    y_diff=function_diff(x_diff)
    #print(df_diff_edge["ZapEnergy"].iloc[-1])
    if options['print'] == True:
        ax3.plot(x_diff,y_diff,color='Green')

    #y_diff_max=np.amax(y_diff,0)
    y_diff_max_index = np.where(y_diff == np.amax(y_diff))
    #print(y_diff_max_index[0])
    edge_shift_diff=float(x_diff[y_diff_max_index])
    print("Edge shift estimated by the differential maximum is "+str(round(edge_shift_diff,5)))
    if options['print'] == True:
        ax3.axvline(x=edge_shift_diff,color="green")
    #print(df_doublediff_edge["ZapEnergy"].iloc[0])
    #ax4.plot(x_doublediff,y_doublediff,color='Green'))


    #fitting double differentiate 
    df_doublediff = df_doublediff[df_doublediff[filenames].notna()]
    d = np.polyfit(df_doublediff_edge["ZapEnergy"],df_doublediff_edge[filenames],2)
    function_doublediff = np.poly1d(d)

    x_doublediff=np.linspace(df_doublediff_edge["ZapEnergy"].iloc[0],df_doublediff_edge["ZapEnergy"].iloc[-1],num=10000)
    y_doublediff=function_doublediff(x_doublediff)

    if options['print'] == True:
        ax4.plot(x_doublediff,y_doublediff,color='Green')

    y_doublediff_zero=find_nearest(y_doublediff,0)
    y_doublediff_zero_index = np.where(y_doublediff == y_doublediff_zero)
    
    edge_shift_doublediff=float(x_doublediff[y_doublediff_zero_index])
   
    print("Edge shift estimated by the double differential zero-point is "+str(round(edge_shift_doublediff,5)))
    if options['print'] == True:
        ax4.axvline(x=edge_shift_doublediff,color="green")
    
