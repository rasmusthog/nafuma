import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import nafuma.auxillary as aux
import nafuma.xanes as xas
import nafuma.xanes.io as io
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
        edge_start = 6.42
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

def post_edge_normalization(path, options={}):

    required_options = ['print']
    default_options = {
        'print': False
    }
    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    
    df_bkgd_sub,filenames,edge = pre_edge_subtraction(path)
    #Defining the end of the pre-edge-region for Mn/Ni, thus start of the edge
    #Implement widget
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

    return df_bkgd_sub, df_postedge