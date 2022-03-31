import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import beamtime.auxillary as aux

def rbkerbest():
    print("ROSENBORG!<3")

#def split_xanes_scan(filename, destination=None):

 #   with open(filename, 'r') as f:


##Better to make a new function that loops through the files, and performing the split_xanes_scan on


def pre_edge_subtraction(df,filenames, options={}):

    required_options = ['edge', 'print']
    default_options = {
        'edge' : 'Mn',
        'print': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    #Defining the end of the pre-edge-region for Mn/Ni, thus start of the edge
    if str(options['edge']) == 'Mn':
        edge_start = 6.45
    if str(options['edge']) == 'Ni':
        edge_start = 8.3


    #making a function to check the difference between values in the list and the defined start of the edge (where background regression will stop):
    absolute_difference_function = lambda list_value : abs(list_value - edge_start)

    #finding the energy data point value that is closest to what I defined as the end of the background
    edge_start_value = min(df["ZapEnergy"], key=absolute_difference_function)

    #Finding what the index of the edge shift end point is
    start_index=df[df["ZapEnergy"]==edge_start_value].index.values[0]

    #Defining x-range for linear background fit, ending at the edge start index
    df_start=df[0:start_index]
    
    #Making a new dataframe, with only the ZapEnergies as the first column
    df_background = pd.DataFrame(df["ZapEnergy"])

    for files in filenames:

    #Fitting linear function to the pre-edge
        d = np.polyfit(df_start["ZapEnergy"],df_start[files],1)
        function_pre = np.poly1d(d)
        
    #making a list, y_pre,so the background will be applied to all ZapEnergy-values
        y_pre=function_pre(df["ZapEnergy"])
        
    #adding a new column in df_background with the y-values of the background
        df_background.insert(1,files,y_pre) 
    
        #Plotting the calculated pre-edge background with the region used for the regression
    
        ###     FOR FIGURING OUT WHERE IT GOES WRONG/WHICH FILES IS CORRUPT
            #ax = df.plot(x = "ZapEnergy",y=files)  
        
    if options['print'] == True:
    #Plotting an example of the edge_start region and the fitted background that will later be subtracted
        ax = df.plot(x = "ZapEnergy",y=filenames[0]) #defining x and y
        plt.axvline(x = edge_start_value) 
        fig = plt.figure(figsize=(15,15))
        df_background.plot(x="ZapEnergy", y=filenames[0],color="Red",ax=ax)
        
    ###################### Subtracting the pre edge from xmap_roi00   ################
    #making a new dataframe to insert the background subtracted intensities
    df_new = pd.DataFrame(df["ZapEnergy"])
    #inserting the pre_edge-background subtracted original xmap_roi00 data

    for files in filenames:
        newintensity_calc=df[files]-df_background[files]
        df_new.insert(1,files,newintensity_calc) 

    if options['print'] == True:
        #Plotting original data (black) and background subtracted data (red)
        ax = df.plot(x = "ZapEnergy",y=filenames[0], color="Black")
        plt.axvline(x = edge_start_value) 
        fig = plt.figure(figsize=(15,15))
        df_new.plot(x="ZapEnergy", y=filenames[0],color="Red",ax=ax)
    return df_new

def post_edge_normalization(df,df_new,filenames, options={}):

    required_options = ['edge', 'print']
    default_options = {
        'edge' : 'Mn',
        'print': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    #Defining the end of the pre-edge-region for Mn/Ni, thus start of the edge
    if str(options['edge']) == 'Mn':
        edge_stop = 6.565
    if str(options['edge']) == 'Ni':
        edge_stop = 8.361

    absolute_difference_function = lambda list_value : abs(list_value - edge_stop) 
    edge_stop_value = min(df_new["ZapEnergy"], key=absolute_difference_function) 
    end_index=df_new[df_new["ZapEnergy"]==edge_stop_value].index.values[0] 
    #Defining x-range for linear fit
    df_fix=df_new
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
        ax = df_new.plot(x = "ZapEnergy",y=filenames) #defining x and y
        plt.axvline(x = edge_stop_value) 
        fig = plt.figure(figsize=(15,15))
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        #print(function_post_list)
        #print(function_post)
        ax = df_new.plot(x = "ZapEnergy",y=filenames, legend=False) #defining x and y
        df_postedge.plot(x="ZapEnergy", y=filenames,color="Green",ax=ax, legend=False)  
        plt.axvline(x = edge_stop_value) 
    
    