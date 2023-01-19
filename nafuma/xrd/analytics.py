import nafuma.auxillary as aux
import nafuma.xrd as xrd
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import datetime
import os

def background_subtracted_peak(data,options):



    default_options = {
    'plot_before_and_after':  False,
    'background_shoulder_left': 0.1, #picking out how far from the peak edge on the left side the background will be extracted from
    'background_shoulder_right': 0.1, #picking out how far from the peak edge on the left side the background will be extracted from
    'peak_interval': False, #Provide an interval [x1,x2]
    'save_dir': 'background_subtracted_peak',
    'background_poly_degree': 1,
    'plot_all_background_fits': False
    }



    diffractogram, wavelength = xrd.io.read_xy(data=data,options=options)   

    if "noheaders" in data['path'][0]:
        filename = os.path.basename(data['path'][0]).split('_noheaders.')[0]
    else:
        filename = os.path.basename(data['path'][0]).split('.')[0]

    options = aux.update_options(options=options, default_options=default_options)
    
    if not options['peak_interval']:
        print("Make sure that you provide an interval defining the 2th-range of the peak")
    
####################################################################################################
#============================ Defining the background in the first dataset ========================
####################################################################################################

    #normalized_background_points_y=[] #one array for each peak, containting the normalized y-values to be multiplied with a factor (the average value of the intensities between the first backround point and the first peak point)
    #background_points_x=[] #one array of all points between the background points boundaries for each peak (one array per peak of interest)
#for j, peak in enumerate(peaks): #for-loop for each of the peaks I want to integrate
    background_left_shoulder_x=[] #for each peak of interest, I hereby fill in the x-values of the background before and after the peak
    background_left_shoulder_y=[] #for each peak of interest, I hereby fill  in the y-values of the background before and after the peak
    background_right_shoulder_x=[] #for each peak of interest, I hereby fill in the x-values of the background before and after the peak
    background_right_shoulder_y=[] #for each peak of interest, I hereby fill  in the y-values of the background before and after the peak
    background_full_x=[]
    data_full_y=[]
    peak_x=[] # for each peak of interest, I hereby fill in the x-values of the peak
    peak_y=[] #for each peak of interest, I hereby fill in the y-values of the peak 
    #background_points_for_normalization=[]
    #x_range=[]

    for i, twotheta in enumerate(diffractogram["2th"]): #using the background start and end points to define the regions of interest
        if options['peak_interval'][0]-options['background_shoulder_left'] < twotheta and twotheta < options['peak_interval'][1]+options['background_shoulder_right']:
            background_full_x.append(twotheta)
            data_full_y.append(diffractogram['I'][i])
            if  twotheta < options['peak_interval'][0]:
                background_left_shoulder_x.append(twotheta)
                background_left_shoulder_y.append(diffractogram["I"][i])
            elif twotheta < options['peak_interval'][1]:
                peak_x.append(twotheta)
                peak_y.append(diffractogram["I"][i])
            elif twotheta < options['peak_interval'][1]+options['background_shoulder_right']:
                background_right_shoulder_x.append(twotheta)
                background_right_shoulder_y.append(diffractogram["I"][i])

    background_shoulders_x=np.concatenate((background_left_shoulder_x, background_right_shoulder_x))
    background_shoulders_y=np.concatenate((background_left_shoulder_y, background_right_shoulder_y))  
    #if backgroundstart[j] < values and values < peakstart[j]: #only picking out the data points between backgroun_start and peak_start, as these will be used to normalize the background
             #   background_points_for_normalization.append(df["intensity_tot"][i])
    #background_points_x.append(background_x)
    #background_normalization_factor = np.average(background_points_for_normalization)
    #normalized_background_y=background_y/background_normalization_factor
    #print('Normalization factor for '+str(peak)+' is '+str(background_normalization_factor))
    #normalized_background_points_y.append(normalized_background_y)
    #print(normalized_background_points_y)
            #x_range.append(values)


            #test


    #print(background_shoulders_x)
    #print(background_shoulders_y)
    d = np.polyfit(background_shoulders_x, background_shoulders_y,options['background_poly_degree'])
    function_background = np.poly1d(d) #Using the values of the background to make a backgroudn (2. deg polynomial)
    
#Applying the fitted function to the twotheta-values of the whole 2-theta region of relevance
    background_y_fitted=function_background(background_full_x)

    #Subtracting the background from the peak to be left with the peak itelf
    data_minus_background = data_full_y - background_y_fitted
    df_peak = pd.DataFrame()
    df_peak['2th']=background_full_x
    df_peak['I_org']=data_full_y
    df_peak['I_BG']=background_y_fitted
    df_peak['I_corr']=data_minus_background

    ################## PLOTTING A CHOSEN NUMBER OF FILES TO DOUBLE CHECK THE VALIDITY OF THE SCRIPT ####################################################################################################################        
    if options['plot_all_background_fits']:

        ax = df_peak.plot(x="2th",y="I_BG")
        
        ax.set_ylim(min(df_peak['I_org'])*0.9,max(df_peak['I_org'])*1.1)
        
        ax.set_xlim(options['peak_interval'][0]-2*(options['background_shoulder_left']+options['background_shoulder_right']),options['peak_interval'][1]+2*(options['background_shoulder_left']+options['background_shoulder_right']))
        ax.set_title(filename)
        diffractogram.plot(x="2th",y="I",ax=ax)
        
        plt.axvline(x = options['peak_interval'][0])
        plt.axvline(x = options['peak_interval'][1])

#######################################################################################################################################


    return diffractogram, df_peak

def find_fwhm_of_peak(x,y,start_values,options):
    #Here the data needs to be two arrays

##############################################################################
# =========== Subtracting the peak from the background
###############################################################################
    default_options={
        'lorentzian':True,
        'voigt': False,
        'pseudovoigt': False,
        #'starting_guess_lor':   [400,   14.1,   0.01],
        #'starting_guess_gauss': [400,   14.1,   0.01],
        'plot_fit': True,
        'lower_bounds': [0,0,0,0,0],
        'lower_bounds_PV': [0,0,0,0,0],
        'upper_bounds': [np.inf,np.inf,np.inf,0],
        'upper_bounds_PV': [np.inf,np.inf,np.inf,1]

    }
    
    options = aux.update_options(options=options, default_options=default_options)

    ################################ Fitting  ###################
    def _1Voigt(x, ampL, center, widL, ampG, sigmaG):
        return (ampG*(1/(sigmaG*(np.sqrt(2*np.pi))))*(np.exp(-((x-center)**2)/((2*sigmaG)**2)))) +\
            ((ampL*widL**2/((x-center)**2+widL**2)) )

    def _1PV(x, I, x0, PV_fwhm, ratio):
        sigma=PV_fwhm/(2*np.sqrt(2*np.log(2)))
        a_G= 1/(sigma*np.sqrt(2*np.pi))
        b_G= 4*np.log(2)/PV_fwhm**2

        GAUSSIAN_PART= a_G*np.exp(-b_G*(x-x0)**2)
        LORENTZIAN_PART= 1/np.pi * (PV_fwhm/2)/((x-x0)**2+(PV_fwhm/2)**2)

        return (I * (ratio * GAUSSIAN_PART+(1-ratio)*LORENTZIAN_PART))

    def _1Lorentzian(x, ampL, center, widL):
        return ((ampL*widL**2/((x-center)**2+widL**2)) )
    
    # defining starting values for the fit
    ampL = start_values[0]
    center = start_values[1]
    widL = start_values[2]

    #producing lots of x-points, to make graph better
    x_fit=np.linspace(min(x),max(x),50*len(x))

    if options['lorentzian']:
        param_bounds=(options['lower_bounds'][:3],options['upper_bounds'][:3])
        popt_lor, pcov_lor = scipy.optimize.curve_fit(_1Lorentzian, x, y, p0=[ampL, center, widL],bounds=param_bounds)
        
        perr_lor = np.sqrt(np.diag(pcov_lor))

        [ampL, center, widL]=popt_lor
        parameters=popt_lor
        errors=perr_lor


        y_fit=  ((ampL*widL**2/((x_fit-center)**2+widL**2)) )
        #     voigt_fit= (ampG1_fitted*(1/(sigmaG1_fitted*(np.sqrt(2*np.pi))))*(np.exp(-((x_values-cenG1_fitted)**2)/((2*sigmaG1_fitted)**2)))) +\
        #   ((ampL1_fitted*widL1_fitted**2/((x_values-cenL1_fitted)**2+widL1_fitted**2)) ) 
    
    
    if options['voigt']:
        ampG = start_values[3]
        #forcing it to be the same peak position, thus no specific center value necessary here
        sigmaG = start_values[-1]
        
        param_bounds=(options['lower_bounds'],options['upper_bounds'])
        popt_voigt, pcov_voigt = scipy.optimize.curve_fit(_1Voigt, x, y, p0=[ampL, center, widL, ampG, sigmaG],bounds=param_bounds)
        
        perr_voigt = np.sqrt(np.diag(pcov_voigt))

        [ampL, center, widL, ampG, sigmaG]=popt_voigt
        parameters=popt_voigt
        errors=perr_voigt


        y_fit= (ampG*(1/(sigmaG*(np.sqrt(2*np.pi))))*(np.exp(-((x_fit-center)**2)/((2*sigmaG)**2)))) +\
            ((ampL*widL**2/((x_fit-center)**2+widL**2)) )

    if options['pseudovoigt']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I= start_values[0]#PeakIntensity 

        x0 = start_values[1]
        PV_fwhm= start_values[2]
        ratio = start_values[3]
                
        param_bounds=(options['lower_bounds_PV'][:4],options['upper_bounds_PV'][:4])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_1PV, x, y, p0=[I,x0,PV_fwhm,ratio],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I,x0,PV_fwhm,ratio]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        sigma=PV_fwhm/(2*np.sqrt(2*np.log(2)))
        a_G= 1/(sigma*np.sqrt(2*np.pi))
        b_G= 4*np.log(2)/PV_fwhm**2

        GAUSSIAN_PART= a_G*np.exp(-b_G*(x_fit-x0)**2)
        LORENTZIAN_PART= 1/np.pi * (PV_fwhm/2)/((x_fit-x0)**2+(PV_fwhm/2)**2)

        y_fit = (I * (ratio * GAUSSIAN_PART+(1-ratio)*LORENTZIAN_PART))

        #y_fit= I * (ratio * 1/(PV_fwhm/(2*np.sqrt(2*np.log(2)))*np.sqrt(2*np.pi))*np.exp(-4*np.log(2)/PV_fwhm**2*(x_fit-x0)**2)+(1-ratio)*1/np.pi * (PV_fwhm/2)/((x_fit-x0)**2+(PV_fwhm/2)**2))
            
    if options['plot_fit']:
        fig = plt.figure(figsize=[40,10])
        ax = plt.axes()
        ax.plot(x,y,'o')
    
    #ax.set_ylim(min(df_peak['I_org'])*0.9,max(df_peak['I_org'])*1.1)
    
    #ax.set_xlim(options['peak_interval'][0]-2*(options['background_shoulder_left']+options['background_shoulder_right']),options['peak_interval'][1]+2*(options['background_shoulder_left']+options['background_shoulder_right']))
    #ax.set_title(filename)
    #diffractogram.plot(x="2th",y="I",ax=ax)
    
    #plt.axvline(x = options['peak_interval'][0])
    #plt.axvline(x = options['peak_interval'][1])
    
    # plt.scatter(x=x,y=y)
        ax.plot(x_fit,y_fit)#,col='b')
    #print(x)
    return parameters, errors

def read_timestamps(data,options):
    default_options={
        'scan_time_s_manual': False, #this is just an example, and should be set to False at default
        'decimals_for_seconds': 3 #3 this is given if data is in milliseconds
       # 'include_microseconds': False
    }

    options = aux.update_options(options=options, default_options=default_options)

    df_time= pd.read_csv(data['path'],delim_whitespace=True, header=0) #relevant headers are "#!Filename", "Date" and "Time"

    #FIXME Start with for-loop to remake all times and dates to one time-date-value that can be easily used to calculate differences (timedelta) 
    # both for scantime-estimate and the relative time axis, which is the main takaway from this function
    #df_time_new=pd.DataFrame(columns=["DateTime"])
    time_stamp_list=[]
    for i,time in enumerate(df_time['Time']):
        time_string=df_time['Time'].iloc[i]
        date_string=df_time['Date'].iloc[i]
        #print(date_string)
        time_stamp=datetime.datetime(
        year        =int(date_string.split('-')[0]),
        month       =int(date_string.split('-')[1]),
        day         =int(date_string.split('-')[2]),
        hour        =int(time_string.split(':')[0]),
        minute      =int(time_string.split(':')[1]),
        second      =int(time_string.split(':')[2].split('.')[0]),
        microsecond =int(time_string.split(':')[2].split('.')[1])*10**options['decimals_for_seconds'])
        
        time_stamp_list.append(time_stamp)
    #print(df_time.columns)
    df_time_new=pd.DataFrame(df_time["#!File_name"].values,columns=["Filename"])
    #df_time_new=pd.DataFrame(time_stamp_list,columns=["DateTime"])
    df_time_new['DateTime'] = time_stamp_list#df_time["#!File_name"].values
    #df_time_new=pd.concat([df_time_new1, df_time_new2], axis=1)
    #print(df_time["#!File_name"])
    
    
    #fixing scan time:
    if options['scan_time_s_manual']:
        scan_time_s=options['scan_time_s_manual']
    else:
        scan_time_s=(df_time_new["DateTime"][20]-df_time_new["DateTime"][0])/20
    
    #FIXME add the print the start time of the experiment
    first_timestamp = df_time_new['DateTime'].iloc[0]
    first_timestamp_corrected = first_timestamp - scan_time_s
    
    print('First data point, corrected for scan time, was from '+str(first_timestamp_corrected))
    relative_time_list_h=[]
    for time in df_time_new["DateTime"]:
        relative_time=time-first_timestamp
        relative_time_in_h=relative_time.days * 24 + relative_time.seconds/3600 + relative_time.microseconds/3600000000
        #print(relative_time_in_h)
        relative_time_list_h.append(relative_time_in_h) 
    
    df_time_new["Relative_time"]=relative_time_list_h
    #print(df_time_new)
    #df_time_new["Filename"]=df_time["#!Filename"]
    return first_timestamp_corrected, df_time_new

    
def find_area_of_peak(x,y,options):
    #x and y should be already background subtracted, using the background-subtracted_peak-function
    default_options={}    
        
    options = aux.update_options(options=options, default_options=default_options)

    #The distance between each x-point is the same throughout the dataset:
    dx=x[1]-x[0]
# Compute the area using the composite trapezoidal rule.
    peak_area = np.trapz(y, dx=dx)

    return peak_area