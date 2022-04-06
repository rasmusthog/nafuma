import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import beamtime.auxillary as aux
import beamtime.xanes as xas
import beamtime.xanes.io as io
def rbkerbest():
    print("ROSENBORG!<3")

#def split_xanes_scan(filename, destination=None):

 #   with open(filename, 'r') as f:


##Better to make a new function that loops through the files, and performing the split_xanes_scan on

#Tryiung to make a function that can decide which edge it is based on the first ZapEnergy-value
def finding_edge(df):
    if 5.9 < df["ZapEnergy"][0] < 6.5:
        edge='Mn'
        return(edge)
    if 8.0 < df["ZapEnergy"][0] < 8.6:
        edge='Ni'
        return(edge)

#def pre_edge_subtraction(df,filenames, options={}):
def test(innmat):
    df_test= xas.io.put_in_dataframe(innmat)
    print(df_test)

def pre_edge_subtraction(path, options={}):
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
    #implement widget
    if edge == 'Mn':
        edge_start = 6.45
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
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))
        df.plot(x = "ZapEnergy",y=filenames[0],ax=ax1) #defining x and y
        plt.axvline(x = max(df_start["ZapEnergy"])) 
        #fig = plt.figure(figsize=(15,15))
        df_bkgd.plot(x="ZapEnergy", y=filenames[0],color="Red",ax=ax1)
        ax1.set_title('Data and fitted background')
    ###################### Subtracting the pre edge from xmap_roi00   ################
    #making a new dataframe to insert the background subtracted intensities
    df_bkgd_sub = pd.DataFrame(df["ZapEnergy"])
    #inserting the pre_edge-background subtracted original xmap_roi00 data

    for files in filenames:
        newintensity_calc=df[files]-df_bkgd[files]
        df_bkgd_sub.insert(1,files,newintensity_calc) 

    if options['print'] == True:
        df.plot(x = "ZapEnergy",y=filenames[0], color="Black", ax=ax2, legend=False)
        plt.axvline(x = max(df_start["ZapEnergy"])) 
        df_bkgd_sub.plot(x="ZapEnergy", y=filenames[0],color="Red",ax=ax2, legend=False)
        ax2.set_title('Data and background-subtracted data')

    return df_bkgd_sub

def post_edge_normalization(df,df_backg_sub,filenames, options={}):

    required_options = ['print']
    default_options = {
        'print': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    
    edge=finding_edge(df)
    #Defining the end of the pre-edge-region for Mn/Ni, thus start of the edge
    #Implement widget
    if edge == 'Mn':
        edge_stop = 6.565
    if edge == 'Ni':
        edge_stop = 8.361

    absolute_difference_function = lambda list_value : abs(list_value - edge_stop) 
    edge_stop_value = min(df_backg_sub["ZapEnergy"], key=absolute_difference_function) 
    end_index=df_backg_sub[df_backg_sub["ZapEnergy"]==edge_stop_value].index.values[0] 
    #Defining x-range for linear fit
    df_fix=df_backg_sub
    df_fix.dropna(inplace=True) #Removing all indexes without any value, as some of the data sets misses the few last data points and fucks up the fit
    df_end=df_fix[end_index:] #The region of interest for the post edge
    #print(df_end)
    #Fitting linear function to the pre-edge using the background corrected intensities to make the post edge fit
    df_postedge = pd.DataFrame(df["ZapEnergy"])

    function_post_list=[]
    for files in filenames: 
        d = np.polyfit(df_end["ZapEnergy"],df_end[files],1)
        function_post = np.poly1d(d)
        y_post=function_post(df["ZapEnergy"])
        function_post_list.append(function_post)
        df_postedge.insert(1,files,y_post) #adding a new column with the y-values of the fitted post edge

    #print(filenames[0])
    #print(df_postedge)    
    #Plotting the background subtracted signal with the post-edge regression line and the start point for the linear regression line
    if options['print'] == True:
        ax = df_backg_sub.plot(x = "ZapEnergy",y=filenames) #defining x and y
        plt.axvline(x = edge_stop_value) 
        fig = plt.figure(figsize=(15,15))
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        #print(function_post_list)
        #print(function_post)
        ax = df_backg_sub.plot(x = "ZapEnergy",y=filenames, legend=False) #defining x and y
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        plt.axvline(x = edge_stop_value) 
 