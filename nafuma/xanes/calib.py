import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import nafuma.auxillary as aux
import nafuma.xanes as xas
import nafuma.xanes.io as io
from scipy.signal import savgol_filter


##Better to make a new function that loops through the files, and performing the split_xanes_scan on

#Trying to make a function that can decide which edge it is based on the first ZapEnergy-value
def finding_edge(df):
    #FIXME add Fe and Co
    if 5.9 < df["ZapEnergy"][0] < 6.5:
        edge='Mn'
        return(edge)
    if 8.0 < df["ZapEnergy"][0] < 8.6:
        edge='Ni'
        return(edge)

def pre_edge_subtraction(path, options={}):
    #FIXME add log-file instead of the troubleshoot-option
    required_options = ['print','troubleshoot']
    default_options = {
        'print': False,
        'troubleshoot': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    filenames = xas.io.get_filenames(path)
    df= xas.io.put_in_dataframe(path)
    edge=finding_edge(df)
    
    #Defining the end of the region used to define the background, thus start of the edge
    
    #######================================================================================================================================================
    #FIXME Trying to implement automatical region determination based on an estimate of the edge shift
    #print(df)
    #estimated_edge_shift, df_diff, df_diff_max = find_pos_maxdiff(df, filenames,options=options)

    #print(estimated_edge_shift)
    #estimated_edge_shift
    ###=========================================================================================================================================================================
    #implement widget
    if edge == 'Mn':
        edge_start = 6.42
        #edge_start = estimated_edge_shift
    if edge == 'Ni':
        edge_start = 8.3

    #making a dataframe only containing the rows that are included in the background subtraction (points lower than where the edge start is defined)
    df_start=df.loc[df["ZapEnergy"] < edge_start]
        
    #Making a new dataframe, with only the ZapEnergies as the first column -> will be filled to include the background data
    df_bkgd = pd.DataFrame(df["ZapEnergy"])

    for files in filenames:

    #Fitting linear function to the background
        d = np.polyfit(df_start["ZapEnergy"],df_start[files],1)
        function_bkgd = np.poly1d(d)
        
    #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        y_bkgd=function_bkgd(df["ZapEnergy"])
        
    #adding a new column in df_background with the y-values of the background
        df_bkgd.insert(1,files,y_bkgd) 
    
        
        if options['troubleshoot'] == True:
        ###     FOR FIGURING OUT WHERE IT GOES WRONG/WHICH FILE IS CORRUPT
            ax = df.plot(x = "ZapEnergy",y=files)  
    #Plotting the calculated pre-edge background with the region used for the regression   
    if options['print'] == True:
    #Plotting an example of the edge_start region and the fitted background that will later be subtracted
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
        df.plot(x="ZapEnergy", y=filenames,color="Black",ax=ax1)
        df_bkgd.plot(x="ZapEnergy", y=filenames,color="Red",ax=ax1)
        plt.axvline(x = max(df_start["ZapEnergy"])) 
        #fig = plt.figure(figsize=(15,15))
        df_bkgd.plot(x="ZapEnergy", y=filenames,color="Red",ax=ax2)
        ax1.set_title('Data and fitted background')
        #Zooming into bacground region to confirm fit and limits looks reasonable
        df.plot(x = "ZapEnergy",y=filenames,ax=ax2) #defining x and y)
        ax2.set_xlim([min(df_start["ZapEnergy"]),max(df_start["ZapEnergy"])+0.01])
        #finding maximum and minimum values in the backgrounds
        min_values=[]
        max_values=[]
        for file in filenames:
            min_values.append(min(df_start[file]))
            max_values.append(max(df_start[file]))
        ax2.set_ylim([min(min_values),max(max_values)])
        plt.axvline(x = max(df_start["ZapEnergy"]))
        #ax2.set_xlim([25, 50])
    ###################### Subtracting the pre edge from xmap_roi00   ################

    #making a new dataframe to insert the background subtracted intensities
    df_bkgd_sub = pd.DataFrame(df["ZapEnergy"])
    #inserting the background subtracted original xmap_roi00 data

    for files in filenames:
        newintensity_calc=df[files]-df_bkgd[files]
        df_bkgd_sub.insert(1,files,newintensity_calc) 

    if options['print'] == True:
        df.plot(x = "ZapEnergy",y=filenames, color="Black", ax=ax3, legend=False)
        #plt.axvline(x = max(df_start["ZapEnergy"])) 
        df_bkgd_sub.plot(x="ZapEnergy", y=filenames,color="Red",ax=ax3, legend=False)
        ax3.set_title('Data and background-subtracted data')

    return df_bkgd_sub,filenames,edge

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
    
