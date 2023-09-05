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
    'background_region': False, #Provide an interval [x1,x2]
    'background_excluded_region': None,
    'peak_interval': False, #Provide an interval [x1,x2]
    'save_dir': 'background_subtracted_peak',
    'background_poly_degree': 1,
    'plot_all_background_fits': False
    }
    options = aux.update_options(options=options, default_options=default_options)
    diffractogram, wavelength = xrd.io.read_xy(data=data,options=options)   

    if "noheaders" in data['path'][0]:
        filename = os.path.basename(data['path'][0]).split('_noheaders.')[0]
    else:
        filename = os.path.basename(data['path'][0]).split('.')[0]

    
    
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
        #if options['peak_interval'][0]-options['background_shoulder_left'] < twotheta and twotheta < options['peak_interval'][1]+options['background_shoulder_right']:
        if options['background_region'][0] < twotheta and twotheta < options['background_region'][1]:
            background_full_x.append(twotheta)
            data_full_y.append(diffractogram['I'][i])
            if  twotheta < options['peak_interval'][0]:
                background_left_shoulder_x.append(twotheta)
                background_left_shoulder_y.append(diffractogram["I"][i])
            elif twotheta < options['peak_interval'][1]:
                peak_x.append(twotheta)
                peak_y.append(diffractogram["I"][i])
            elif twotheta < options['background_region'][1]:
                background_right_shoulder_x.append(twotheta)
                background_right_shoulder_y.append(diffractogram["I"][i])

    background_shoulders_x=np.concatenate((background_left_shoulder_x, background_right_shoulder_x))
    background_shoulders_y=np.concatenate((background_left_shoulder_y, background_right_shoulder_y))  
    #print(options['background_excluded_region'])
    
    if options['background_excluded_region']:
        for region in options['background_excluded_region']:
            #print(region)
            for i, xval in enumerate(background_shoulders_x):
                if region[0] < xval < region[1]:
                    #because these are numpy arrays I am going via normal arrays to remove the values within the excluded regions              
                    background_x_to_list=background_shoulders_x.tolist()
                    background_y_to_list=background_shoulders_y.tolist()
                    #finding the index of the x-value in the background-array
                    index_xval=background_x_to_list.index(xval)
                    #removing the y-value with the same index
                    background_y_to_list.pop(index_xval)
                    #going back to a numpy array again
                    background_shoulders_y = np.asarray(background_y_to_list)

                    #then removing the xval fro the x-values
                    background_x_to_list.remove(xval)
                    background_shoulders_x = np.asarray(background_x_to_list)

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
    #df_peak['log_BG']=np.log(df_peak['I_BG'])
    
    ################## PLOTTING A CHOSEN NUMBER OF FILES TO DOUBLE CHECK THE VALIDITY OF THE SCRIPT ####################################################################################################################        
    if options['plot_all_background_fits']:
        fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(30,4))
        #fig.suptitle("Background removed data with fit") # or plt.suptitle('Main title')
        #axes[0].plot(x,y,'o')
        #axes[0].plot(x_fit,y_fit)
        df_peak.plot(x="2th",y="I_BG",ax=axes[0])
        
        axes[0].set_ylim(min(df_peak['I_org'])*0.9,min(df_peak['I_org'])*1.5)
        
        #ax.set_xlim(options['peak_interval'][0]-2*(options['background_shoulder_left']+options['background_shoulder_right']),options['peak_interval'][1]+2*(options['background_shoulder_left']+options['background_shoulder_right']))
        axes[0].set_xlim(options['background_region'][0]*0.99,options['background_region'][1]*1.01)
        axes[0].set_title(filename)
        diffractogram.plot(x="2th",y="I",ax=axes[0])
        
        axes[0].axvline(x = options['peak_interval'][0],c='r')
        axes[0].axvline(x = options['peak_interval'][1],c='r')
        axes[0].axvline(x=options['background_region'][0],c='m')
        axes[0].axvline(x=options['background_region'][1],c='m')
        #A plot of only the background subtraction for the ordered peak, to makes sure that is okay from visual inspection
     #   df_peak.plot(x="2th",y="I_BG",ax=axes[1])
     #   axes[1].set_title(filename)

    #    axes[1].set_xlim(options['background_region'][0]*0.99,options['background_region'][0]+(options['background_region'][1]-options['background_region'][0])*0.3)

        #[[35.18,35,76],[36.34,36.59],[36.61,36.79]]
        if options['background_excluded_region']:
            for i in range(len(options['background_excluded_region'])):
                plt.axvline(x = options['background_excluded_region'][i][0],c='g')
                plt.axvline(x = options['background_excluded_region'][i][1],c='g')
        #plt.axvline(x = options['background_excluded_region'][0])
        #plt.axvline(x = options['background_excluded_region'][0])


#######################################################################################################################################


    return diffractogram, df_peak

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

def _2PV(x, I_1, x0_1, PV_fwhm_1, ratio_1, I_2, x0_2, PV_fwhm_2, ratio_2):
    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2))

def _5PV(x, I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5):

    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
    a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
    b_G_3= 4*np.log(2)/PV_fwhm_3**2

    GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x-x0_3)**2)
    LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x-x0_3)**2+(PV_fwhm_3/2)**2)
    
    sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
    a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
    b_G_4= 4*np.log(2)/PV_fwhm_4**2

    GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x-x0_4)**2)
    LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x-x0_4)**2+(PV_fwhm_4/2)**2)

    sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
    a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
    b_G_5= 4*np.log(2)/PV_fwhm_5**2

    GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x-x0_5)**2)
    LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x-x0_5)**2+(PV_fwhm_5/2)**2)

    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5))

def _6PV(x, I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6):

    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
    a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
    b_G_3= 4*np.log(2)/PV_fwhm_3**2

    GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x-x0_3)**2)
    LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x-x0_3)**2+(PV_fwhm_3/2)**2)
    
    sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
    a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
    b_G_4= 4*np.log(2)/PV_fwhm_4**2

    GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x-x0_4)**2)
    LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x-x0_4)**2+(PV_fwhm_4/2)**2)

    sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
    a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
    b_G_5= 4*np.log(2)/PV_fwhm_5**2

    GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x-x0_5)**2)
    LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x-x0_5)**2+(PV_fwhm_5/2)**2)

    sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
    a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
    b_G_6= 4*np.log(2)/PV_fwhm_6**2

    GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x-x0_6)**2)
    LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x-x0_6)**2+(PV_fwhm_6/2)**2)
    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))

def _7PV(x, I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6,I_7,x0_7,PV_fwhm_7,ratio_7):

    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
    a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
    b_G_3= 4*np.log(2)/PV_fwhm_3**2

    GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x-x0_3)**2)
    LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x-x0_3)**2+(PV_fwhm_3/2)**2)
    
    sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
    a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
    b_G_4= 4*np.log(2)/PV_fwhm_4**2

    GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x-x0_4)**2)
    LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x-x0_4)**2+(PV_fwhm_4/2)**2)

    sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
    a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
    b_G_5= 4*np.log(2)/PV_fwhm_5**2

    GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x-x0_5)**2)
    LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x-x0_5)**2+(PV_fwhm_5/2)**2)

    sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
    a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
    b_G_6= 4*np.log(2)/PV_fwhm_6**2

    GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x-x0_6)**2)
    LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x-x0_6)**2+(PV_fwhm_6/2)**2)

    sigma_7=PV_fwhm_7/(2*np.sqrt(2*np.log(2)))
    a_G_7= 1/(sigma_7*np.sqrt(2*np.pi))
    b_G_7= 4*np.log(2)/PV_fwhm_7**2

    GAUSSIAN_PART_7= a_G_7*np.exp(-b_G_7*(x-x0_7)**2)
    LORENTZIAN_PART_7= 1/np.pi * (PV_fwhm_7/2)/((x-x0_7)**2+(PV_fwhm_7/2)**2)

    x0_6 < x0_7

    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))

def _5PV_constraints(x,    
                  I_1,  x0_1,   PV_fwhm_1,  ratio_1,
                  I_2,          PV_fwhm_2,  ratio_2,
                  I_3,  x0_3,   PV_fwhm_3,  ratio_3,
                  I_5,  x0_5,   PV_fwhm_5,  ratio_5,
                  I_6,          PV_fwhm_6,  ratio_6):

    #Testing to force the same phase to have same lattice parameter
    #Next step: making sure I_3/I_4 = I_6/I_7 an
    x0_2=_peak_pos_calc(x0_5,[2,2,2],[3,1,1])
    #x0_4=_peak_pos_calc(x0_1,[3,1,0],[3,1,1])
    x0_6=_peak_pos_calc(x0_3,[3,1,1],[2,2,2])
    #x0_7=_peak_pos_calc(x0_1,[3,1,0],[2,2,2])
    #I_7 = I_6 * I_4/I_3

    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
    a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
    b_G_3= 4*np.log(2)/PV_fwhm_3**2

    GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x-x0_3)**2)
    LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x-x0_3)**2+(PV_fwhm_3/2)**2)
    
    sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
    a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
    b_G_5= 4*np.log(2)/PV_fwhm_5**2

    GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x-x0_5)**2)
    LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x-x0_5)**2+(PV_fwhm_5/2)**2)

    sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
    a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
    b_G_6= 4*np.log(2)/PV_fwhm_6**2

    GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x-x0_6)**2)
    LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x-x0_6)**2+(PV_fwhm_6/2)**2)

    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))


def _7PV_constraints(x,    
                  I_1,  x0_1,   PV_fwhm_1,  ratio_1,
                  I_2,          PV_fwhm_2,  ratio_2,
                  I_3,  x0_3,   PV_fwhm_3,  ratio_3,
                  I_4,          PV_fwhm_4,  ratio_4,
                  I_5,  x0_5,   PV_fwhm_5,  ratio_5,
                  I_6,          PV_fwhm_6,  ratio_6,
                                PV_fwhm_7,  ratio_7):

    #Testing to force the same phase to have same lattice parameter
    #Next step: making sure I_3/I_4 = I_6/I_7 an
    x0_2=_peak_pos_calc(x0_5,[2,2,2],[3,1,1])
    x0_4=_peak_pos_calc(x0_1,[3,1,0],[3,1,1])
    x0_6=_peak_pos_calc(x0_3,[3,1,1],[2,2,2])
    x0_7=_peak_pos_calc(x0_1,[3,1,0],[2,2,2])
    I_7 = I_6 * I_4/I_3

    sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
    a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
    b_G_1= 4*np.log(2)/PV_fwhm_1**2

    GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x-x0_1)**2)
    LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x-x0_1)**2+(PV_fwhm_1/2)**2)

    sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
    a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
    b_G_2= 4*np.log(2)/PV_fwhm_2**2

    GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x-x0_2)**2)
    LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x-x0_2)**2+(PV_fwhm_2/2)**2)

    sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
    a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
    b_G_3= 4*np.log(2)/PV_fwhm_3**2

    GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x-x0_3)**2)
    LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x-x0_3)**2+(PV_fwhm_3/2)**2)
    
    sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
    a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
    b_G_4= 4*np.log(2)/PV_fwhm_4**2

    GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x-x0_4)**2)
    LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x-x0_4)**2+(PV_fwhm_4/2)**2)

    sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
    a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
    b_G_5= 4*np.log(2)/PV_fwhm_5**2

    GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x-x0_5)**2)
    LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x-x0_5)**2+(PV_fwhm_5/2)**2)

    sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
    a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
    b_G_6= 4*np.log(2)/PV_fwhm_6**2

    GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x-x0_6)**2)
    LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x-x0_6)**2+(PV_fwhm_6/2)**2)

    sigma_7=PV_fwhm_7/(2*np.sqrt(2*np.log(2)))
    a_G_7= 1/(sigma_7*np.sqrt(2*np.pi))
    b_G_7= 4*np.log(2)/PV_fwhm_7**2

    GAUSSIAN_PART_7= a_G_7*np.exp(-b_G_7*(x-x0_7)**2)
    LORENTZIAN_PART_7= 1/np.pi * (PV_fwhm_7/2)/((x-x0_7)**2+(PV_fwhm_7/2)**2)

    return (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))


def _1Lorentzian(x, ampL, center, widL):
    return ((ampL*widL**2/((x-center)**2+widL**2)) )

def _1Gaussian(x, I, x0, fwhm_G):
    return (I * np.exp(-(x - x0) ** 2 / (2 * fwhm_G ** 2)))

def find_fwhm_of_peak(x,y,start_values,options):
    #Here the data needs to be two arrays

##############################################################################
# =========== Subtracting the peak from the background
###############################################################################
    default_options={
        'lorentzian':False,
        'voigt': False,
        'pseudovoigt': False,
        'doublePV': False,
        'gaussian': False,
        'cluster_PV': False,
        'cluster_without_RS1': False,
        'cluster_without_RS2': False,
        'cluster_without_RS': False,
        'cluster_split_PV': False,
        'cluster_fullsplit_PV': False,
        'PV_cluster_split': False,
        'PV_cluster_nosplit':False,
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

    if options['gaussian']:
        param_bounds=(options['lower_bounds'][:3],options['upper_bounds'][:3])
        popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_1Gaussian, x, y, p0=[start_values[0], start_values[1], start_values[2]],bounds=param_bounds)
        
        perr_gauss = np.sqrt(np.diag(pcov_gauss))

        [I, x0, fwhm_G]=popt_gauss
        parameters=popt_gauss
        errors=perr_gauss


        y_fit=  (I * np.exp(-(x_fit - x0) ** 2 / (2 * fwhm_G ** 2)))#((ampL*widL**2/((x_fit-center)**2+widL**2)) )
    
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
        I_1= start_values[0]#PeakIntensity 

        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]
                
        param_bounds=(options['lower_bounds_PV'][:4],options['upper_bounds_PV'][:4])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_1PV, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1]=popt_PV
        parameters=popt_PV
        errors=perr_PV
        
        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1))

        #y_fit= I * (ratio * 1/(PV_fwhm/(2*np.sqrt(2*np.log(2)))*np.sqrt(2*np.pi))*np.exp(-4*np.log(2)/PV_fwhm**2*(x_fit-x0)**2)+(1-ratio)*1/np.pi * (PV_fwhm/2)/((x_fit-x0)**2+(PV_fwhm/2)**2))

    if options['doublePV']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]
    
        I_2= start_values[4]#PeakIntensity 
        x0_2 = start_values[5]
        PV_fwhm_2= start_values[6]
        ratio_2 = start_values[7]
        #print('lower bounds below')
        #print(options['lower_bounds_PV'][:4])
        lower_bound=options['lower_bounds_PV'][:4]
        lower_bounds=lower_bound
        
        for i in range(len(lower_bound)):
            lower_bounds.append(lower_bound[i])

        upper_bound=options['upper_bounds_PV'][:4]
        upper_bounds=upper_bound
        for i in range(len(upper_bound)):
            upper_bounds.append(upper_bound[i])
        #print("lower bounds is "+str(lower_bounds))
        #print("upper bounds is "+str(upper_bounds))
        #lower_bounds=np.concatenate(options['lower_bounds_PV'][:4],options['lower_bounds_PV'][:4])
        #upper_bounds=np.concatenate(options['upper_bounds_PV'][:4],options['upper_bounds_PV'][:4])
        param_bounds=(lower_bounds,upper_bounds)
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_2PV, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2))

    if options['cluster_PV']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]

        I_2= start_values[4]#PeakIntensity 
        x0_2 = start_values[5]
        PV_fwhm_2= start_values[6]
        ratio_2 = start_values[7]

        I_3= start_values[8]#PeakIntensity 
        x0_3 = start_values[9]
        PV_fwhm_3= start_values[10]
        ratio_3 = start_values[11]
        
        I_4= start_values[12]#PeakIntensity 
        x0_4 = start_values[13]
        PV_fwhm_4= start_values[14]
        ratio_4 = start_values[15]

        I_5= start_values[16]#PeakIntensity 
        x0_5 = start_values[17]
        PV_fwhm_5= start_values[18]
        ratio_5 = start_values[19]
        
    #######============== defining the lower bound of the fit, deciding it should be defined separately for all the five peaks rather than having the same bounds
   #
   #      lower_bound=options['lower_bounds_PV'][:4]
   #     lower_bounds=lower_bound
        
 #       for i in range(len(lower_bound)):
  #          lower_bounds.append(lower_bound[i])

#        upper_bound=options['upper_bounds_PV'][:4]
#        upper_bounds=upper_bound
#        for i in range(len(upper_bound)):
#            upper_bounds.append(upper_bound[i])

        param_bounds=(options['lower_bounds_PV'],options['upper_bounds_PV'])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_5PV, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
        a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
        b_G_3= 4*np.log(2)/PV_fwhm_3**2

        GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x_fit-x0_3)**2)
        LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x_fit-x0_3)**2+(PV_fwhm_3/2)**2)
        
        sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
        a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
        b_G_4= 4*np.log(2)/PV_fwhm_4**2

        GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x_fit-x0_4)**2)
        LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x_fit-x0_4)**2+(PV_fwhm_4/2)**2)

        sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
        a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
        b_G_5= 4*np.log(2)/PV_fwhm_5**2

        GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x_fit-x0_5)**2)
        LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x_fit-x0_5)**2+(PV_fwhm_5/2)**2)

        if options['cluster_without_RS1']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5))
        elif options['cluster_without_RS2']:
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1)+ I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5))
        elif options['cluster_without_RS']:
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5))  
        else:
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5))
        
    if options['cluster_split_PV']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]

        I_2= start_values[4]#PeakIntensity 
        x0_2 = start_values[5]
        PV_fwhm_2= start_values[6]
        ratio_2 = start_values[7]

        I_3= start_values[8]#PeakIntensity 
        x0_3 = start_values[9]
        PV_fwhm_3= start_values[10]
        ratio_3 = start_values[11]
        
        I_4= start_values[12]#PeakIntensity 
        x0_4 = start_values[13]
        PV_fwhm_4= start_values[14]
        ratio_4 = start_values[15]

        I_5= start_values[16]#PeakIntensity 
        x0_5 = start_values[17]
        PV_fwhm_5= start_values[18]
        ratio_5 = start_values[19]
        
        I_6= start_values[20]#PeakIntensity 
        x0_6 = start_values[21]
        PV_fwhm_6= start_values[22]
        ratio_6 = start_values[23]
    #######============== defining the lower bound of the fit, deciding it should be defined separately for all the five peaks rather than having the same bounds
   #
   #      lower_bound=options['lower_bounds_PV'][:4]
   #     lower_bounds=lower_bound
        
 #       for i in range(len(lower_bound)):
  #          lower_bounds.append(lower_bound[i])

#        upper_bound=options['upper_bounds_PV'][:4]
#        upper_bounds=upper_bound
#        for i in range(len(upper_bound)):
#            upper_bounds.append(upper_bound[i])

        param_bounds=(options['lower_bounds_PV'],options['upper_bounds_PV'])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_6PV, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
        a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
        b_G_3= 4*np.log(2)/PV_fwhm_3**2

        GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x_fit-x0_3)**2)
        LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x_fit-x0_3)**2+(PV_fwhm_3/2)**2)
        
        sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
        a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
        b_G_4= 4*np.log(2)/PV_fwhm_4**2

        GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x_fit-x0_4)**2)
        LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x_fit-x0_4)**2+(PV_fwhm_4/2)**2)

        sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
        a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
        b_G_5= 4*np.log(2)/PV_fwhm_5**2

        GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x_fit-x0_5)**2)
        LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x_fit-x0_5)**2+(PV_fwhm_5/2)**2)

        sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
        a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
        b_G_6= 4*np.log(2)/PV_fwhm_6**2

        GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x_fit-x0_6)**2)
        LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x_fit-x0_6)**2+(PV_fwhm_6/2)**2)

        if options['cluster_without_RS1']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))
        else: #this means if cluster should have RS1 and peak with split
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))

    if options['cluster_fullsplit_PV']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]

        I_2= start_values[4]#PeakIntensity 
        x0_2 = start_values[5]
        #x0_2=_222_to_311(x0_5)
        PV_fwhm_2= start_values[6]
        ratio_2 = start_values[7]

        I_3= start_values[8]#PeakIntensity 
        x0_3 = start_values[9]
        PV_fwhm_3= start_values[10]
        ratio_3 = start_values[11]
        
        I_4= start_values[12]#PeakIntensity 
        x0_4 = start_values[13]
        PV_fwhm_4= start_values[14]
        ratio_4 = start_values[15]

        I_5= start_values[16]#PeakIntensity 
        x0_5 = start_values[17]
        PV_fwhm_5= start_values[18]
        ratio_5 = start_values[19]
        
        I_6= start_values[20]#PeakIntensity 
        x0_6 = start_values[21]
        PV_fwhm_6= start_values[22]
        ratio_6 = start_values[23]

        I_7= start_values[24]#PeakIntensity 
        x0_7 = start_values[25]
        PV_fwhm_7= start_values[26]
        ratio_7 = start_values[27]
    #######============== defining the lower bound of the fit, deciding it should be defined separately for all the five peaks rather than having the same bounds
   #
   #      lower_bound=options['lower_bounds_PV'][:4]
   #     lower_bounds=lower_bound
        
 #       for i in range(len(lower_bound)):
  #          lower_bounds.append(lower_bound[i])

#        upper_bound=options['upper_bounds_PV'][:4]
#        upper_bounds=upper_bound
#        for i in range(len(upper_bound)):
#            upper_bounds.append(upper_bound[i])

        param_bounds=(options['lower_bounds_PV'],options['upper_bounds_PV'])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_7PV, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6,I_7,x0_7,PV_fwhm_7,ratio_7],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,x0_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,x0_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,x0_6,PV_fwhm_6,ratio_6,I_7,x0_7,PV_fwhm_7,ratio_7]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
        a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
        b_G_3= 4*np.log(2)/PV_fwhm_3**2

        GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x_fit-x0_3)**2)
        LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x_fit-x0_3)**2+(PV_fwhm_3/2)**2)
        
        sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
        a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
        b_G_4= 4*np.log(2)/PV_fwhm_4**2

        GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x_fit-x0_4)**2)
        LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x_fit-x0_4)**2+(PV_fwhm_4/2)**2)

        sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
        a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
        b_G_5= 4*np.log(2)/PV_fwhm_5**2

        GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x_fit-x0_5)**2)
        LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x_fit-x0_5)**2+(PV_fwhm_5/2)**2)

        sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
        a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
        b_G_6= 4*np.log(2)/PV_fwhm_6**2

        GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x_fit-x0_6)**2)
        LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x_fit-x0_6)**2+(PV_fwhm_6/2)**2)

        sigma_7=PV_fwhm_7/(2*np.sqrt(2*np.log(2)))
        a_G_7= 1/(sigma_7*np.sqrt(2*np.pi))
        b_G_7= 4*np.log(2)/PV_fwhm_7**2

        GAUSSIAN_PART_7= a_G_7*np.exp(-b_G_7*(x_fit-x0_7)**2)
        LORENTZIAN_PART_7= 1/np.pi * (PV_fwhm_7/2)/((x_fit-x0_7)**2+(PV_fwhm_7/2)**2)

        if options['cluster_without_RS1']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
        elif options['cluster_without_RS2']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1)+ I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
        elif options['cluster_without_RS']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
      
        else: #this means if cluster should have RS1 and peak with split
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))

    if options['PV_cluster_split']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]

        I_2= start_values[4]#PeakIntensity 
        #x0_2 = start_values[5]
        PV_fwhm_2= start_values[5]
        ratio_2 = start_values[6]

        I_3= start_values[7]#PeakIntensity 
        x0_3 = start_values[8]
        PV_fwhm_3= start_values[9]
        ratio_3 = start_values[10]
        
        I_4= start_values[11]#PeakIntensity 
        #x0_4 = start_values[13]
        PV_fwhm_4= start_values[12]
        ratio_4 = start_values[13]

        I_5= start_values[14]#PeakIntensity 
        x0_5 = start_values[15]
        PV_fwhm_5= start_values[16]
        ratio_5 = start_values[17]
        
        I_6= start_values[18]#PeakIntensity 
        #x0_6 = start_values[21]
        PV_fwhm_6= start_values[19]
        ratio_6 = start_values[20]

        #I_7= I_6*I_4/I_3 #PeakIntensity 
        #x0_7 = start_values[25]
        PV_fwhm_7= start_values[21]
        ratio_7 = start_values[22]
    #######============== defining the lower bound of the fit, deciding it should be defined separately for all the five peaks rather than having the same bounds
   #
   #      lower_bound=options['lower_bounds_PV'][:4]
   #     lower_bounds=lower_bound
        
 #       for i in range(len(lower_bound)):
  #          lower_bounds.append(lower_bound[i])

#        upper_bound=options['upper_bounds_PV'][:4]
#        upper_bounds=upper_bound
#        for i in range(len(upper_bound)):
#            upper_bounds.append(upper_bound[i])

        param_bounds=(options['lower_bounds_PV'],options['upper_bounds_PV'])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_7PV_constraints, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,PV_fwhm_6,ratio_6,PV_fwhm_7,ratio_7],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_4,PV_fwhm_4,ratio_4,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,PV_fwhm_6,ratio_6,PV_fwhm_7,ratio_7]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        x0_2=_peak_pos_calc(x0_5,[2,2,2],[3,1,1])
        x0_4=_peak_pos_calc(x0_1,[3,1,0],[3,1,1])
        x0_6=_peak_pos_calc(x0_3,[3,1,1],[2,2,2])
        x0_7=_peak_pos_calc(x0_1,[3,1,0],[2,2,2])
        I_7 = I_6 * I_4/I_3

        parameters=np.insert(parameters,21,x0_7)
        parameters=np.insert(parameters,21,I_7)
        parameters=np.insert(parameters,19,x0_6)
        parameters=np.insert(parameters,12,x0_4)
        parameters=np.insert(parameters,5,x0_2)
    
        x0_1_err=errors[1]
        x0_3_err=errors[8]
        x0_5_err=errors[15]
        I_4_err = errors[11]

        x0_2_err=x0_5_err
        x0_4_err=x0_1_err
        x0_6_err=x0_3_err
        x0_7_err=x0_1_err
        #estimated error for the I_7 intensity, just scaling with the I_4 intensity error:
        I_7_err = I_4_err * I_7/I_4

        errors=np.insert(errors,21,x0_7_err)
        errors=np.insert(errors,21,I_7_err)
        errors=np.insert(errors,19,x0_6_err)
        errors=np.insert(errors,12,x0_4_err)
        errors=np.insert(errors,5,x0_2_err)
#######################==========================================

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
        a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
        b_G_3= 4*np.log(2)/PV_fwhm_3**2

        GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x_fit-x0_3)**2)
        LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x_fit-x0_3)**2+(PV_fwhm_3/2)**2)
        
        sigma_4=PV_fwhm_4/(2*np.sqrt(2*np.log(2)))
        a_G_4= 1/(sigma_4*np.sqrt(2*np.pi))
        b_G_4= 4*np.log(2)/PV_fwhm_4**2

        GAUSSIAN_PART_4= a_G_4*np.exp(-b_G_4*(x_fit-x0_4)**2)
        LORENTZIAN_PART_4= 1/np.pi * (PV_fwhm_4/2)/((x_fit-x0_4)**2+(PV_fwhm_4/2)**2)

        sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
        a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
        b_G_5= 4*np.log(2)/PV_fwhm_5**2

        GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x_fit-x0_5)**2)
        LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x_fit-x0_5)**2+(PV_fwhm_5/2)**2)

        sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
        a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
        b_G_6= 4*np.log(2)/PV_fwhm_6**2

        GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x_fit-x0_6)**2)
        LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x_fit-x0_6)**2+(PV_fwhm_6/2)**2)

        sigma_7=PV_fwhm_7/(2*np.sqrt(2*np.log(2)))
        a_G_7= 1/(sigma_7*np.sqrt(2*np.pi))
        b_G_7= 4*np.log(2)/PV_fwhm_7**2

        GAUSSIAN_PART_7= a_G_7*np.exp(-b_G_7*(x_fit-x0_7)**2)
        LORENTZIAN_PART_7= 1/np.pi * (PV_fwhm_7/2)/((x_fit-x0_7)**2+(PV_fwhm_7/2)**2)

        if options['cluster_without_RS1']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
        elif options['cluster_without_RS2']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1)+ I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
        elif options['cluster_without_RS']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))
      
        else: #this means if cluster should have RS1 and peak with split
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_4 * (ratio_4 * GAUSSIAN_PART_4+(1-ratio_4)*LORENTZIAN_PART_4)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6)+ I_7 * (ratio_7 * GAUSSIAN_PART_7+(1-ratio_7)*LORENTZIAN_PART_7))

    if options['PV_cluster_nosplit']:
        #https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
        I_1= start_values[0]#PeakIntensity 
        x0_1 = start_values[1]
        PV_fwhm_1= start_values[2]
        ratio_1 = start_values[3]

        I_2= start_values[4]#PeakIntensity 
        #x0_2 = start_values[5]
        PV_fwhm_2= start_values[5]
        ratio_2 = start_values[6]

        I_3= start_values[7]#PeakIntensity 
        x0_3 = start_values[8]
        PV_fwhm_3= start_values[9]
        ratio_3 = start_values[10]

        I_5= start_values[11]#PeakIntensity 
        x0_5 = start_values[12]
        PV_fwhm_5= start_values[13]
        ratio_5 = start_values[14]
        
        I_6= start_values[15]#PeakIntensity 
        #x0_6 = start_values[21]
        PV_fwhm_6= start_values[16]
        ratio_6 = start_values[17]

    #######============== defining the lower bound of the fit, deciding it should be defined separately for all the five peaks rather than having the same bounds
   #
   #      lower_bound=options['lower_bounds_PV'][:4]
   #     lower_bounds=lower_bound
        
 #       for i in range(len(lower_bound)):
  #          lower_bounds.append(lower_bound[i])

#        upper_bound=options['upper_bounds_PV'][:4]
#        upper_bounds=upper_bound
#        for i in range(len(upper_bound)):
#            upper_bounds.append(upper_bound[i])

        param_bounds=(options['lower_bounds_PV'],options['upper_bounds_PV'])
        popt_PV, pcov_PV = scipy.optimize.curve_fit(_5PV_constraints, x, y, p0=[I_1,x0_1,PV_fwhm_1,ratio_1,I_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,PV_fwhm_6,ratio_6],bounds=param_bounds)
        
        perr_PV = np.sqrt(np.diag(pcov_PV))

        [I_1,x0_1,PV_fwhm_1,ratio_1,I_2,PV_fwhm_2,ratio_2,I_3,x0_3,PV_fwhm_3,ratio_3,I_5,x0_5,PV_fwhm_5,ratio_5,I_6,PV_fwhm_6,ratio_6]=popt_PV
        parameters=popt_PV
        errors=perr_PV

        x0_2=_peak_pos_calc(x0_5,[2,2,2],[3,1,1])
        #x0_4=_peak_pos_calc(x0_1,[3,1,0],[3,1,1])
        x0_6=_peak_pos_calc(x0_3,[3,1,1],[2,2,2])
        #x0_7=_peak_pos_calc(x0_1,[3,1,0],[2,2,2])
        #I_7 = I_6 * I_4/I_3

        #parameters=np.insert(parameters,21,x0_7)
        #parameters=np.insert(parameters,21,I_7)
        parameters=np.insert(parameters,16,x0_6)
        #parameters=np.insert(parameters,12,x0_4)
        parameters=np.insert(parameters,5,x0_2)
        x0_3_err=errors[8]
        x0_5_err=errors[12]
        #I_4_err = errors[11]

        x0_2_err=x0_5_err
        #x0_4_err=x0_1_err
        x0_6_err=x0_3_err
        #x0_7_err=x0_1_err
        #estimated error for the I_7 intensity, just scaling with the I_4 intensity error:
        #I_7_err = I_4_err * I_7/I_4

        #errors=np.insert(errors,21,x0_7_err)
        #errors=np.insert(errors,21,I_7_err)
        errors=np.insert(errors,16,x0_6_err)
        #errors=np.insert(errors,12,x0_4_err)
        errors=np.insert(errors,5,x0_2_err)
#######################==========================================

        sigma_1=PV_fwhm_1/(2*np.sqrt(2*np.log(2)))
        a_G_1= 1/(sigma_1*np.sqrt(2*np.pi))
        b_G_1= 4*np.log(2)/PV_fwhm_1**2

        GAUSSIAN_PART_1= a_G_1*np.exp(-b_G_1*(x_fit-x0_1)**2)
        LORENTZIAN_PART_1= 1/np.pi * (PV_fwhm_1/2)/((x_fit-x0_1)**2+(PV_fwhm_1/2)**2)

        sigma_2=PV_fwhm_2/(2*np.sqrt(2*np.log(2)))
        a_G_2= 1/(sigma_2*np.sqrt(2*np.pi))
        b_G_2= 4*np.log(2)/PV_fwhm_2**2

        GAUSSIAN_PART_2= a_G_2*np.exp(-b_G_2*(x_fit-x0_2)**2)
        LORENTZIAN_PART_2= 1/np.pi * (PV_fwhm_2/2)/((x_fit-x0_2)**2+(PV_fwhm_2/2)**2)

        sigma_3=PV_fwhm_3/(2*np.sqrt(2*np.log(2)))
        a_G_3= 1/(sigma_3*np.sqrt(2*np.pi))
        b_G_3= 4*np.log(2)/PV_fwhm_3**2

        GAUSSIAN_PART_3= a_G_3*np.exp(-b_G_3*(x_fit-x0_3)**2)
        LORENTZIAN_PART_3= 1/np.pi * (PV_fwhm_3/2)/((x_fit-x0_3)**2+(PV_fwhm_3/2)**2)
        
        sigma_5=PV_fwhm_5/(2*np.sqrt(2*np.log(2)))
        a_G_5= 1/(sigma_5*np.sqrt(2*np.pi))
        b_G_5= 4*np.log(2)/PV_fwhm_5**2

        GAUSSIAN_PART_5= a_G_5*np.exp(-b_G_5*(x_fit-x0_5)**2)
        LORENTZIAN_PART_5= 1/np.pi * (PV_fwhm_5/2)/((x_fit-x0_5)**2+(PV_fwhm_5/2)**2)

        sigma_6=PV_fwhm_6/(2*np.sqrt(2*np.log(2)))
        a_G_6= 1/(sigma_6*np.sqrt(2*np.pi))
        b_G_6= 4*np.log(2)/PV_fwhm_6**2

        GAUSSIAN_PART_6= a_G_6*np.exp(-b_G_6*(x_fit-x0_6)**2)
        LORENTZIAN_PART_6= 1/np.pi * (PV_fwhm_6/2)/((x_fit-x0_6)**2+(PV_fwhm_6/2)**2)
        
        if options['cluster_without_RS1']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))
        elif options['cluster_without_RS2']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1)+ I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))
        elif options['cluster_without_RS']:
            #removing any contribution from I2, the RS1-peak
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))
      
        else: #this means if cluster should have RS1 and peak with split
            y_fit = (I_1 * (ratio_1 * GAUSSIAN_PART_1+(1-ratio_1)*LORENTZIAN_PART_1) + I_2 * (ratio_2 * GAUSSIAN_PART_2+(1-ratio_2)*LORENTZIAN_PART_2)+ I_3 * (ratio_3 * GAUSSIAN_PART_3+(1-ratio_3)*LORENTZIAN_PART_3)+ I_5 * (ratio_5 * GAUSSIAN_PART_5+(1-ratio_5)*LORENTZIAN_PART_5)+ I_6 * (ratio_6 * GAUSSIAN_PART_6+(1-ratio_6)*LORENTZIAN_PART_6))



    if options['plot_fit']:
        fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(30,4))
        #fig.suptitle("Background removed data with fit") # or plt.suptitle('Main title')
        axes[0].plot(x,y,'o')
        axes[0].plot(x_fit,y_fit)

        axes[1].plot(x,np.cbrt(y),'o')
        axes[1].plot(x_fit,np.cbrt(y_fit))
        
        axes[2].plot(x,y,'o')
        axes[2].plot(x_fit,y_fit)
        ##TEST


        axes[0].axvline(parameters[1]) #center position
        x_fwhm1=[parameters[1]-parameters[2]/2,parameters[1]+parameters[2]/2] #plotting the width of the fwhm
        #FIXME generalize for other functions than PV:
#        sigma=PV_fwhm/(2*np.sqrt(2*np.log(2)))
 #       a_G= 1/(sigma*np.sqrt(2*np.pi))      
        #y=[(I*a_G)/2,(I*a_G)/2]
        
        #y_fwhm=[max(y_fit)/2,max(y_fit)/2]

        #is this useful? yes for PV it should be useful
        y1_max=I_1 * (ratio_1 * a_G_1+(1-ratio_1)*(1/(np.pi*PV_fwhm_1/2)))
        y_fwhm1=[y1_max/2,y1_max/2]
        axes[0].plot(x_fwhm1,y_fwhm1,c='g')
        axes[1].plot(x_fwhm1,np.cbrt(y_fwhm1),c='g')
        
        #Trying to plot the ordering-peak
        #axes[2].plot(x_fit,y1_max)
        axes[2].set_ylim(0,y1_max*1.5)

        if options['doublePV']:
            plt.axvline(parameters[5]) #center position
            x_fwhm2=[parameters[5]-parameters[6]/2,parameters[5]+parameters[6]/2] #plotting the width of the fwhm
            #FIXME generalize for other functions than PV:
    #        sigma=PV_fwhm/(2*np.sqrt(2*np.log(2)))
    #       a_G= 1/(sigma*np.sqrt(2*np.pi))      
            #y=[(I*a_G)/2,(I*a_G)/2]

            y2_max=I_2 * (ratio_2 * a_G_2+(1-ratio_2)*(1/(np.pi*PV_fwhm_2/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm2=[y2_max/2,y2_max/2]
            axes[0].plot(x_fwhm2,y_fwhm2,c='g')
            ##TEST
        if options['cluster_PV'] or options['PV_cluster_nosplit']:
            if not options['cluster_without_RS1'] or not options['cluster_without_RS']:    
                plt.axvline(parameters[5]) #center position
                x_fwhm2=[parameters[5]-parameters[6]/2,parameters[5]+parameters[6]/2] #plotting the width of the fwhm
                y2_max=I_2 * (ratio_2 * a_G_2+(1-ratio_2)*(1/(np.pi*PV_fwhm_2/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm2=[y2_max/2,y2_max/2]
                axes[0].plot(x_fwhm2,y_fwhm2,c='b')
                axes[1].plot(x_fwhm2,np.cbrt(y_fwhm2),c='b')

            plt.axvline(parameters[9]) #center position
            x_fwhm3=[parameters[9]-parameters[10]/2,parameters[9]+parameters[10]/2] #plotting the width of the fwhm
            y3_max=I_3 * (ratio_3 * a_G_3+(1-ratio_3)*(1/(np.pi*PV_fwhm_3/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm3=[y3_max/2,y3_max/2]
            axes[0].plot(x_fwhm3,y_fwhm3,c='r')
            axes[1].plot(x_fwhm3,np.cbrt(y_fwhm3),c='r')

            if not options['cluster_without_RS2'] or not options['cluster_without_RS']: 
                plt.axvline(parameters[13]) #center position
                x_fwhm5=[parameters[13]-parameters[14]/2,parameters[13]+parameters[14]/2] #plotting the width of the fwhm
                y5_max=I_5 * (ratio_5 * a_G_5+(1-ratio_5)*(1/(np.pi*PV_fwhm_5/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm5=[y5_max/2,y5_max/2]
                axes[0].plot(x_fwhm5,y_fwhm5,c='b')
                axes[1].plot(x_fwhm5,np.cbrt(y_fwhm5),c='b')

            plt.axvline(parameters[17]) #center position
            x_fwhm6=[parameters[17]-parameters[18]/2,parameters[17]+parameters[18]/2] #plotting the width of the fwhm
            y6_max=I_6 * (ratio_6 * a_G_6+(1-ratio_6)*(1/(np.pi*PV_fwhm_6/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm6=[y6_max/2,y6_max/2]
            axes[0].plot(x_fwhm6,y_fwhm6,c='r')
            axes[1].plot(x_fwhm6,np.cbrt(y_fwhm6),c='r')
        '''
        if options['cluster_split_PV']:
            plt.axvline(parameters[5]) #center position
            x_fwhm2=[parameters[5]-parameters[6]/2,parameters[5]+parameters[6]/2] #plotting the width of the fwhm
            y2_max=I_2 * (ratio_2 * a_G_2+(1-ratio_2)*(1/(np.pi*PV_fwhm_2/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm2=[y2_max/2,y2_max/2]
            axes[0].plot(x_fwhm2,y_fwhm2,c='g')
            axes[1].plot(x_fwhm2,np.cbrt(y_fwhm2),c='g')

            if not options['cluster_without_RS1'] or not options['cluster_without_RS']:
                plt.axvline(parameters[9]) #center position
                x_fwhm3=[parameters[9]-parameters[10]/2,parameters[9]+parameters[10]/2] #plotting the width of the fwhm
                y3_max=I_3 * (ratio_3 * a_G_3+(1-ratio_3)*(1/(np.pi*PV_fwhm_3/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm3=[y3_max/2,y3_max/2]
                axes[0].plot(x_fwhm3,y_fwhm3,c='b')
                axes[1].plot(x_fwhm3,np.cbrt(y_fwhm3),c='b')

            plt.axvline(parameters[13]) #center position
            x_fwhm4=[parameters[13]-parameters[14]/2,parameters[13]+parameters[14]/2] #plotting the width of the fwhm
            y4_max=I_4 * (ratio_4 * a_G_4+(1-ratio_4)*(1/(np.pi*PV_fwhm_4/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm4=[y4_max/2,y4_max/2]
            axes[0].plot(x_fwhm4,y_fwhm4,c='g')
            axes[1].plot(x_fwhm4,np.cbrt(y_fwhm4),c='g')

            if not options['cluster_without_RS2'] or not options['cluster_without_RS']:
                plt.axvline(parameters[17]) #center position
                x_fwhm5=[parameters[17]-parameters[18]/2,parameters[17]+parameters[18]/2] #plotting the width of the fwhm
                y5_max=I_5 * (ratio_5 * a_G_5+(1-ratio_5)*(1/(np.pi*PV_fwhm_5/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm5=[y5_max/2,y5_max/2]
                axes[0].plot(x_fwhm5,y_fwhm5,c='b')
                axes[1].plot(x_fwhm5,np.cbrt(y_fwhm5),c='b')

            plt.axvline(parameters[21]) #center position
            x_fwhm6=[parameters[21]-parameters[22]/2,parameters[21]+parameters[22]/2] #plotting the width of the fwhm
            y6_max=I_6 * (ratio_6 * a_G_6+(1-ratio_6)*(1/(np.pi*PV_fwhm_6/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm6=[y6_max/2,y6_max/2]
            axes[0].plot(x_fwhm6,y_fwhm6,c='g')
            axes[1].plot(x_fwhm6,np.cbrt(y_fwhm6),c='g')
        '''    
            
        if options['cluster_fullsplit_PV'] or options['PV_cluster_split']:
            if not options['cluster_without_RS1'] or not options['cluster_without_RS']:
                plt.axvline(parameters[5]) #center position
                x_fwhm2=[parameters[5]-parameters[6]/2,parameters[5]+parameters[6]/2] #plotting the width of the fwhm
                y2_max=I_2 * (ratio_2 * a_G_2+(1-ratio_2)*(1/(np.pi*PV_fwhm_2/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm2=[y2_max/2,y2_max/2]
                axes[0].plot(x_fwhm2,y_fwhm2,c='b')
                axes[1].plot(x_fwhm2,np.cbrt(y_fwhm2),c='b')

            plt.axvline(parameters[9]) #center position
            x_fwhm3=[parameters[9]-parameters[10]/2,parameters[9]+parameters[10]/2] #plotting the width of the fwhm
            y3_max=I_3 * (ratio_3 * a_G_3+(1-ratio_3)*(1/(np.pi*PV_fwhm_3/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm3=[y3_max/2,y3_max/2]
            axes[0].plot(x_fwhm3,y_fwhm3,c='r')
            axes[1].plot(x_fwhm3,np.cbrt(y_fwhm3),c='r')

            plt.axvline(parameters[13]) #center position
            x_fwhm4=[parameters[13]-parameters[14]/2,parameters[13]+parameters[14]/2] #plotting the width of the fwhm
            y4_max=I_4 * (ratio_4 * a_G_4+(1-ratio_4)*(1/(np.pi*PV_fwhm_4/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm4=[y4_max/2,y4_max/2]
            axes[0].plot(x_fwhm4,y_fwhm4,c='g')
            axes[1].plot(x_fwhm4,np.cbrt(y_fwhm4),c='g')

            if not options['cluster_without_RS2'] or not options['cluster_without_RS']:
                plt.axvline(parameters[17]) #center position
                x_fwhm5=[parameters[17]-parameters[18]/2,parameters[17]+parameters[18]/2] #plotting the width of the fwhm
                y5_max=I_5 * (ratio_5 * a_G_5+(1-ratio_5)*(1/(np.pi*PV_fwhm_5/2))) #trying to reduce the second peak at x=x0_2
                y_fwhm5=[y5_max/2,y5_max/2]
                axes[0].plot(x_fwhm5,y_fwhm5,c='b')
                axes[1].plot(x_fwhm5,np.cbrt(y_fwhm5),c='b')

            plt.axvline(parameters[21]) #center position
            x_fwhm6=[parameters[21]-parameters[22]/2,parameters[21]+parameters[22]/2] #plotting the width of the fwhm
            y6_max=I_6 * (ratio_6 * a_G_6+(1-ratio_6)*(1/(np.pi*PV_fwhm_6/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm6=[y6_max/2,y6_max/2]
            axes[0].plot(x_fwhm6,y_fwhm6,c='r')
            axes[1].plot(x_fwhm6,np.cbrt(y_fwhm6),c='r')

            plt.axvline(parameters[25]) #center position
            x_fwhm7=[parameters[25]-parameters[26]/2,parameters[25]+parameters[26]/2] #plotting the width of the fwhm
            y7_max=I_7 * (ratio_7 * a_G_7+(1-ratio_7)*(1/(np.pi*PV_fwhm_7/2))) #trying to reduce the second peak at x=x0_2
            y_fwhm7=[y7_max/2,y7_max/2]
            axes[0].plot(x_fwhm7,y_fwhm7,c='g')
            axes[1].plot(x_fwhm7,np.cbrt(y_fwhm7),c='g')
        
        
            ##TEST
            #ax.set_ylim(min(df_peak['I_org'])*0.9,max(df_peak['I_org'])*1.1)
    
    #ax.set_xlim(options['peak_interval'][0]-2*(options['background_shoulder_left']+options['background_shoulder_right']),options['peak_interval'][1]+2*(options['background_shoulder_left']+options['background_shoulder_right']))
    #ax.set_title(filename)
    #diffractogram.plot(x="2th",y="I",ax=ax)
    
    #plt.axvline(x = options['peak_interval'][0])
    #plt.axvline(x = options['peak_interval'][1])
    
    # plt.scatter(x=x,y=y)

        
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
        time_stamp=datetime.datetime(
        year        =int(date_string.split('-')[0]),
        month       =int(date_string.split('-')[1]),
        day         =int(date_string.split('-')[2]),
        hour        =int(time_string.split(':')[0]),
        minute      =int(time_string.split(':')[1]),
        second      =int(time_string.split(':')[2].split('.')[0]),
        microsecond =int(time_string.split(':')[2].split('.')[1])*10**options['decimals_for_seconds'])
        
        time_stamp_list.append(time_stamp)
    df_time_new=pd.DataFrame(df_time["#!File_name"].values,columns=["Filename"])
    #df_time_new=pd.DataFrame(time_stamp_list,columns=["DateTime"])
    df_time_new['DateTime'] = time_stamp_list#df_time["#!File_name"].values
    #df_time_new=pd.concat([df_time_new1, df_time_new2], axis=1)
    
    
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
        relative_time_list_h.append(relative_time_in_h) 
    
    df_time_new["Relative_time"]=relative_time_list_h

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

def _DC1(x, ad, bd, cd):
        
        return (np.sqrt(ad * np.cos(x*np.pi/180)**4 + bd * np.cos(x*np.pi/180)**2 + cd))

def find_fit_parameters_for_peak(chosen_peaks,options): #can add wavelength if it will be used on data from another beam time at some point
    default_options={
        'temperature': False, #FIXME this is not implemented, but would be a great addition if running on variable-temperature data-sets, using some kind of thermal expansion constant
        'wavelength': False 
                   }
    options = aux.update_options(options=options, default_options=default_options)
    
    start_values_list=[]
    #start_values_list2=[] #for regions of overlapping peaks, such as disorder peak which is splitting over time. Ideally only starting values is necessary, and the background used and the peak region could be the same.
    peak_range_list=[]
    #peak_range_list2=[] #this might or might not be necessary, depending on how the fitting goes.
    background_range_list=[]
    #number_of_excluded_regions_list=[]
    BG_poly_degree_list=[]
    excluded_background_range_list=[] #FIXME add excluded background regions in the analyze-background_subtracted-function
    #df_peaks=pd.DataFrame(columns=chosen_peaks)
    if "ord1" in chosen_peaks:
        peak_center= 14.095
        start_values_list.append([0.1, peak_center, 0.2434226, 0.5])
        peak_range_list.append(         [13.9,14.3])
        background_range_list.append(   [13.5,14.35])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values

    if "Pt3" in chosen_peaks:
        peak_center=42.431
        start_values_list.append([100, peak_center, 0.1, 0.8]) #bad first guesses, fix this
        peak_range_list.append([42.3,42.56])
        background_range_list.append([42.25,42.75])
        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #three excluded regions in the background region
        #df_peaks["Pt3"]=start_values
    #print(excluded_background_range_list)
        BG_poly_degree_list.append(2)
    if "char_split" in chosen_peaks: #char as in a characteristic peak that both disordered and ordered phase will have
        #FIXME Have to add the possibility to have 
        peak_center=36.10 #43.35
        peak_center2=36.17 #43.45
        peak_end=36.27
        peak_range_list.append([35.82, peak_end])
        background_range_list.append([35.72,37.25])
        start_values_list.append([34.63923911, peak_center,  0.0847148,   0.92987947, 34.63923911, peak_center,  0.0847148,   0.92987947]) 
        #start_values_list2.append([34.63923911, peak_center2,  0.0847148,   0.92987947])
        #number_of_excluded_regions_list.append(3) #three excluded regions in the background region
        excluded_background_range_list.append([[peak_end,37.18]])#[[35.23,35.76],[36.34,37.22]]) #include this for certain peaks that has peaks in close proximity and influences the background
        BG_poly_degree_list.append(1)
    if "char" in chosen_peaks: #char as in a characteristic peak that both disordered and ordered phase will have
        #FIXME Have to add the possibility to have 
        peak_center=36.16 #43.35
        peak_end=36.27
        peak_range_list.append([35.82, peak_end])
        background_range_list.append([35.72,37.25])#[35.76,36.87])
        start_values_list.append([34.63923911, peak_center,  0.0847148,   0.92987947])#9.61189651 42.43286743  0.07566345  0.816174 
        #start_values_list2.append(None)
        #number_of_excluded_regions_list.append(3) #three excluded regions in the background region
        excluded_background_range_list.append([[peak_end,37.18]])#[[35.18,35.74],[36.34,36.59],[36.61,36.79]])#([,,) #include this for certain peaks that has peaks in close proximity and influences the background
        BG_poly_degree_list.append(1)
    if options['temperature']:
        original_wavelength=0.6390512 #in Å
        #multiply all values with a temperature dependent factor and change values accordingly
        print("tempereture dependence not yet implemented")
    return BG_poly_degree_list,start_values_list,background_range_list, peak_range_list, excluded_background_range_list


def WL_translate(twotheta_original,wavelength_original,wavelength_new):
    #new_peak= 2*np.arcsin(reference_wavelength/data["wavelength"][0]) * 180/np.pi
    #only including peaks that are actually in the data set
    twotheta_new=2*np.arcsin((wavelength_new/wavelength_original)*np.sin(twotheta_original/2 * np.pi/180))*180/np.pi
    return twotheta_new
'''
def _222_to_311(peak_pos_222):
    #input is the peak position (2th) of the (2 2 2) rock salt peak, and output is estimated pos of the (3 1 1).
    #calculated using braggs law and the relation between a and d for cubic structures (a=1/sqrt(h^2+k^2+l^2))
    peak_pos_311=2*180/np.pi*np.arcsin(np.sqrt(11/12)*np.sin(peak_pos_222/2*np.pi/180))
    return peak_pos_311
'''
def _peak_pos_calc(pos_0,miller_0,miller_new):
    #input is the peak position (2th) of the (2 2 2) rock salt peak, and output is estimated pos of the (3 1 1).
    #calculated using braggs law and the relation between a and d for cubic structures (a=1/sqrt(h^2+k^2+l^2))
    d_ratio=np.sqrt(miller_new[0]**2+miller_new[1]**2+miller_new[2]**2)/np.sqrt(miller_0[0]**2+miller_0[1]**2+miller_0[2]**2)
    peak_pos_new=2*180/np.pi*np.arcsin(d_ratio*np.sin(pos_0/2*np.pi/180))
    return peak_pos_new


def _lattice_param_calc(pos,pos_err,miller,wavelength):
    #Braggs law
    theta=pos/2 * np.pi/180 
    theta_err = pos_err/2 * np.pi/180 
    d = wavelength/(2*np.sin(theta))
    #d_err = wavelength/(2*np.cos(theta_err))
    #lattice params and miller indices for cubic systems
    a = d * np.sqrt(miller[0]**2+miller[1]**2+miller[2]**2)
    a_err = np.sqrt(wavelength/(2)     *     (miller[0]**2+miller[1]**2+miller[2]**2)      *    (-np.cos(theta)/np.sin(theta)**2)**2)   * theta_err 
    return a,a_err

def find_fit_parameters_for_peak_general(chosen_peaks,options): #can add wavelength if it will be used on data from another beam time at some point
    #FIXME add peak splitting of the cluster, 
    default_options={
        'temperature': False, #FIXME this is not implemented, but would be a great addition if running on variable-temperature data-sets, using some kind of thermal expansion constant
        'wavelength': False 
                   }
    
    options = aux.update_options(options=options, default_options=default_options)
    WL=0.6390512
    T=25
    if options['wavelength']:
        WL_new=options['wavelength']
    else:
        print("No wavelength given in options!")

    if options['temperature']:
        T_diff=options['temperature']-T
    else:
        T_diff = 0

    start_values_list=[]
    #start_values_list2=[] #for re_genergions of overlapping peaks, such as disorder peak which is splitting over time. Ideally only starting values is necessary, and the background used and the peak region could be the same.
    peak_range_list=[]
    #peak_range_list2=[] #this might or might not be necessary, depending on how the fitting goes.
    background_range_list=[]
    #number_of_excluded_regions_list=[]
    BG_poly_degree_list=[]
    excluded_background_range_list=[] #FIXME add excluded background regions in the analyze-background_subtracted-function
    #df_peaks=pd.DataFrame(columns=chosen_peaks)

    lower_bounds_list=[]
    upper_bounds_list=[]
    twoth_exp=0.00001 * T_diff #approximately the increase in twotheta per degree K above RT (just an approx, without any basis from theory)
    
    slack=0.12
    slack_RS=0.05
    slack_ord=0.015

    if "ord1" in chosen_peaks:
        peak_center= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        start_values_list.append([0.1, peak_center, 0.2434226, 0.5])
        lower_bounds_list.append([0,       peak_center-slack_ord,    0.01,  0.5])
        upper_bounds_list.append([20,       peak_center+slack_ord,    0.3,    1])
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(14.4*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.5*(1-twoth_exp),WL,WL_new),WL_translate(14.45*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values

    ##Trying to implement a cluster of peaks, have to consider whether sub1 and sub2 must be separated into two peaks as well to acount for peak splitting
    if "cluster" in chosen_peaks:
        peak_center_ord= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        peak_center_RS1= WL_translate(14.67*(1-twoth_exp),WL,WL_new)#14.59, 14.67
        peak_center_sub1= WL_translate(14.9*(1-twoth_exp),WL,WL_new)
        peak_center_RS2 = WL_translate(15.32*(1-twoth_exp),WL,WL_new)
        peak_center_sub2= WL_translate(15.53*(1-twoth_exp),WL,WL_new)
        start_values_list.append(
            [0,    peak_center_ord,    0.15,  1,
             0,    peak_center_RS1,    0.05,  1, 
             30,   peak_center_sub1,   0.05,  1,
             0,  peak_center_RS2,    0.05,  1,
             10,   peak_center_sub2,   0.05,  1]
        )

        lower_bounds_list.append(
            [0,       peak_center_ord-slack_ord,    0.01,  0.5,
             0,       peak_center_RS1-slack_RS,    0.01,  0.5, 
             0,       peak_center_sub1-slack,   0.01,  0.5,
             0,       peak_center_RS2-slack_RS,    0.01,  0.5,
             0,       peak_center_sub2-slack,   0.01,  0.5]
        )   
        upper_bounds_list.append(
            [20,       peak_center_ord+slack_ord,    0.3,    1,
             20,       peak_center_RS1+slack_RS,    0.2,    1, 
             1200,       peak_center_sub1+slack,   0.3,  1,
             150,       peak_center_RS2+slack_RS,    0.3,    1,
             600,       peak_center_sub2+slack,   0.3, 1]
        )   
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(15.8*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.5*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values 
    #slack_sub=0.05
    '''
    if "cluster_split" in chosen_peaks:
        peak_center_ord= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        peak_center_RS1= WL_translate(14.67*(1-twoth_exp),WL,WL_new)#14.59, 14.67
        peak_center_sub1= WL_translate(14.9*(1-twoth_exp),WL,WL_new)
        peak_center_RS2 = WL_translate(15.32*(1-twoth_exp),WL,WL_new)
        peak_center_sub2_dis= WL_translate(15.52*(1-twoth_exp),WL,WL_new)
        peak_center_sub2_ord= WL_translate(15.56*(1-twoth_exp),WL,WL_new)

        start_values_list.append(
            [0,    peak_center_ord,    0.15,  1,
             0,    peak_center_RS1,    0.05,  1, 
             30,   peak_center_sub1,   0.05,  1,
             0,  peak_center_RS2,    0.05,  1,
            10, peak_center_sub2_dis, 0.05,1,
             0,   peak_center_sub2_ord,   0.05,  1]
        )

        lower_bounds_list.append(
            [0,       peak_center_ord-slack_ord,    0.01,  0.5,
             0,       peak_center_RS1-slack_RS,    0.01,  0.5, 
             0,       peak_center_sub1-slack,   0.01,  0.5,
             0,       peak_center_RS2-slack_RS,    0.01,  0.5,
             0,       peak_center_sub2_dis-slack,   0.01,  0.5,
             0,       peak_center_sub2_ord-slack_ord,   0.01,  0.5]
        )   
        upper_bounds_list.append(
            [20,       peak_center_ord+slack_ord,    0.3,    1,
             20,       peak_center_RS1+slack_RS,    0.2,    1, 
             1200,       peak_center_sub1+slack,   0.3,  1,
             150,       peak_center_RS2+slack_RS,    0.3,    1,
             600,       peak_center_sub2_dis+slack,   0.3, 1,
             500,       peak_center_sub2_ord+slack_ord,   0.3, 1]
        )   
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(15.8*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.6*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values 
    '''
    if "cluster_fullsplit" in chosen_peaks:
        peak_center_ord= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        peak_center_RS1= WL_translate(14.65*(1-twoth_exp),WL,WL_new)#14.59, 14.67
        peak_center_sub1_dis= WL_translate(14.88*(1-twoth_exp),WL,WL_new)
        peak_center_sub1_ord= WL_translate(14.92*(1-twoth_exp),WL,WL_new)
        peak_center_RS2 = WL_translate(15.32*(1-twoth_exp),WL,WL_new)
        peak_center_sub2_dis= WL_translate(15.52*(1-twoth_exp),WL,WL_new)
        peak_center_sub2_ord= WL_translate(15.56*(1-twoth_exp),WL,WL_new)

        start_values_list.append(
            [0,    peak_center_ord,    0.15,  1,
             0,    peak_center_RS1,    0.05,  1, 
             10,   peak_center_sub1_dis,   0.05,  1,
             0,   peak_center_sub1_ord,   0.05,  1,
             0,  peak_center_RS2,    0.05,  1,
            3, peak_center_sub2_dis, 0.05,1,
             0,   peak_center_sub2_ord,   0.05,  1]
        )

        lower_bounds_list.append(
            [0,       peak_center_ord-slack_ord,    0.05,  0,
             0,       peak_center_RS1-slack_RS,    0.05,  0, 
             0,       peak_center_sub1_dis-slack,   0.02,  0,
             0,       peak_center_sub1_ord-slack_ord,   0.02,  0,
             0,       peak_center_RS2-slack_RS,    0.05,  0,
             0,       peak_center_sub2_dis-slack,   0.02,  0,
             0,       peak_center_sub2_ord-slack_ord,   0.02,  0]
        )   
        upper_bounds_list.append(
            [20,       peak_center_ord+slack_ord,    0.3,    1,
             20,       peak_center_RS1+slack_RS,    0.2,    1, 
             1200,       peak_center_sub1_dis+slack,   0.18,  1,
             1000,       peak_center_sub1_ord+slack_ord,   0.18, 1,
             150,       peak_center_RS2+slack_RS,    0.25,    1,
             600,       peak_center_sub2_dis+slack,   0.18, 1,
             500,       peak_center_sub2_ord+slack_ord,   0.18, 1]
        )   
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(15.8*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.6*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values 
    
    if 'PV_cluster_split' in chosen_peaks:
        peak_center_ord= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        
        peak_center_sub1_dis= WL_translate(14.80*(1-twoth_exp),WL,WL_new) #14.88
        #peak_center_sub1_ord= _peak_pos_calc(peak_center_ord,[3,1,0],[3,1,1])# WL_translate(14.92*(1-twoth_exp),WL,WL_new)
        peak_center_RS2 = WL_translate(15.32*(1-twoth_exp),WL,WL_new)
        #peak_center_sub2_dis= _peak_pos_calc(peak_center_sub1_dis,[3,1,1],[2,2,2])#WL_translate(15.52*(1-twoth_exp),WL,WL_new)
        #peak_center_sub2_ord= _peak_pos_calc(peak_center_ord,[3,1,0],[2,2,2])#WL_translate(15.56*(1-twoth_exp),WL,WL_new)

        #peak_center_RS1= _peak_pos_calc(peak_center_RS2,[2,2,2],[3,1,1])#WL_translate(14.65*(1-twoth_exp),WL,WL_new)#14.59, 14.67
        
        start_values_list.append(
            [   0,    peak_center_ord,         0.15,   1,
                0,                             0.05,   1, 
                5,   peak_center_sub1_dis,    0.05,   1,
                0,                             0.05,   1,
                0,  peak_center_RS2,           0.05,   1,
                2,                             0.05,   1,
                                               0.05,   1]
        )

        lower_bounds_list.append(
            [0,       peak_center_ord-slack_ord,    0.05,  0,
             0,                                     0.05,  0, 
             0,       peak_center_sub1_dis-slack,   0.02,  0.8,
             0,                                     0.02,  0.8,
             0,       peak_center_RS2-slack_RS,    0.05,  0,
             0,                                     0.02,  0.8,
                                                    0.02,  0.8]
        )   
        upper_bounds_list.append(
            [20,       peak_center_ord+slack_ord,    0.3,    1,
             20,                                    0.2,    1, 
             1200,       peak_center_sub1_dis+slack,   0.18,  1,
             1000,                                      0.18, 1,
             150,       peak_center_RS2+slack_RS,    0.25,    1,
             600,                                   0.18, 1,
                                                    0.18, 1]
        )   
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(15.8*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.6*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values 
    
    if 'PV_cluster_nosplit' in chosen_peaks:
        peak_center_ord= WL_translate(14.2*(1-twoth_exp),WL,WL_new)
        peak_center_sub1= WL_translate(14.88*(1-twoth_exp),WL,WL_new)
        peak_center_RS2 = WL_translate(15.32*(1-twoth_exp),WL,WL_new)
        
        start_values_list.append(
            [   0,    peak_center_ord,         0.15,   1,
                0,                             0.05,   1, 
                10,   peak_center_sub1,         0.05,   1,
                0,    peak_center_RS2,           0.05,   1,
                3,                             0.05,   1]
        )

        lower_bounds_list.append(
            [0,       peak_center_ord-slack_ord,    0.05,  0,
             0,                                     0.05,  0, 
             0,       peak_center_sub1-slack,   0.02,  0.75,
             0,       peak_center_RS2-slack_RS,    0.05,  0,
             0,                                     0.02,  0.75,]
        )   
        upper_bounds_list.append(
            [20,       peak_center_ord+slack_ord,    0.3,    1,
             20,                                    0.2,    1, 
             1200,       peak_center_sub1+slack,    0.18,  1,
             150,       peak_center_RS2+slack_RS,    0.25,    1,
             600,                                   0.18, 1,]
        )   
        peak_range_list.append(         [WL_translate(13.9*(1-twoth_exp),WL,WL_new),WL_translate(15.805*(1-twoth_exp),WL,WL_new)])
        #background_range_list.append(   [WL_translate(13.6*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        background_range_list.append(   [WL_translate(13.5*(1-twoth_exp),WL,WL_new),WL_translate(15.95*(1-twoth_exp),WL,WL_new)])
        BG_poly_degree_list.append(2)        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]])#13.1,13.5]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #no excluded regions
        #df_peaks["ord1"]=start_values 

    if "Pt3" in chosen_peaks:
        peak_center=WL_translate(42.431,WL,WL_new)
        start_values_list.append([100, peak_center, 0.1, 0.8]) #bad first guesses, fix this
        peak_range_list.append([WL_translate(42.3,WL,WL_new),WL_translate(42.56,WL,WL_new)])
        background_range_list.append([WL_translate(42.25,WL,WL_new),WL_translate(42.75,WL,WL_new)])
        #start_values_list2.append(None)
        excluded_background_range_list.append([[0,0]]) #include this for certain peaks that has peaks in close proximity
        #number_of_excluded_regions_list.append(0) #three excluded regions in the background region
        #df_peaks["Pt3"]=start_values
    #print(excluded_background_range_list)
        BG_poly_degree_list.append(2)
    if "char_split" in chosen_peaks: #char as in a characteristic peak that both disordered and ordered phase will have
        #FIXME Have to add the possibility to have 
        peak_center=WL_translate(36.10,WL,WL_new) #43.35
        peak_center2=WL_translate(36.17,WL,WL_new) #43.45
        peak_end=WL_translate(36.27,WL,WL_new)
        peak_range_list.append([WL_translate(35.82,WL,WL_new), peak_end])
        background_range_list.append([WL_translate(35.72,WL,WL_new),WL_translate(37.25,WL,WL_new)])
        start_values_list.append([34.63923911, peak_center,  0.0847148,   0.92987947, 34.63923911, peak_center2,  0.0847148,   0.92987947]) 
        #start_values_list2.append([34.63923911, peak_center2,  0.0847148,   0.92987947])
        #number_of_excluded_regions_list.append(3) #three excluded regions in the background region
        excluded_background_range_list.append([[peak_end,WL_translate(37.18,WL,WL_new)]])#[[35.23,35.76],[36.34,37.22]]) #include this for certain peaks that has peaks in close proximity and influences the background
        BG_poly_degree_list.append(1)

    return lower_bounds_list,upper_bounds_list, BG_poly_degree_list,start_values_list,background_range_list, peak_range_list, excluded_background_range_list

def finding_instrumental_peak_broadening(data,options):
    #This function might not be necessary, as I can use the DC1-function and refine data with TOPAS to find proper parameters. But parts of this can be used for other things.
    #FIXME
    #use one WL as the standard and plot in all the peak positions
    #translate these peak positions to the wavelength given
    #only keep the peak positions that are less than the biggest x-value in the file (so not too many peaks are included and we get an error)
    #use this to find the start_values-list for each of the peaks, by using som standard distance on each side of the peak and for the background. This should be rather simple for lab6.
    #Find the fwhm of all peaks, and have an option to print the fwhm as a function of angle
    #Fit these to a function (fwhm vs 2theta), to extract the estimated fwhm on the specific angle of interest
    #As a return value the fwhm of the relevant angle is given as e.g. fwhm_instrumental
    #Then as this value is given to the function working on the actual data, any fwhm value equal to or below this should give an error as that should be impossible and something is obviously wrong
    default_options = {
       # 'lorentzian': True,
        'peak_center': 0,
        'plot_instrumental_broadening': True,
        'plot_every_nth': 10
        #'background_poly_degree': 1
    }
    
    options = aux.update_options(options=options, default_options=default_options)
    
    diffractogram, wavelength = xrd.io.read_xy(data=data,options=options)  

    start_values_list=[]
    peak_range_list=[]
    background_range_list=[]
    
    reference_peak_centers=[8.81678204487931, 12.484732202720007, 15.305588264031849, 17.687216986202827, 19.793091426749534, 21.703371161909057, 25.113360089664365,
    26.663096159547234, 28.13383083445108, 29.537041574463057, 30.88185501694746, 32.178100946275805, 33.42858705018042, 35.81271399776499, 36.953180739661214, 
    38.06362325243035, 39.14737970910566, 40.20672810099413, 41.24399752773122, 42.26064266026674, 44.24359987984366, 45.20783351961699, 46.15589302740396,
    47.086112710436716, 48.908189795887864, 49.80081533788387, 51.55182963558574, 52.4] #all values except the last is found through peak fitting
    reference_wavelength=0.6390512
    peak_centers=[]
    #peak_centers=reference_peak_centers


    #value to be used to pick out peak-region and background region
    peak_shoulders=0.4 
    background_shoulders=0.45
    
#FIXME translation of wavelength 
    #translating to correct wavelength
    for peak in reference_peak_centers:
        #new_peak= 2*np.arcsin(reference_wavelength/data["wavelength"][0]) * 180/np.pi
        #only including peaks that are actually in the data set
        new_peak=2*np.arcsin((data["wavelength"][0]/reference_wavelength)*np.sin(peak/2 * np.pi/180))*180/np.pi
        #print(new_peak)
        if new_peak > min(diffractogram["2th"])+background_shoulders and new_peak < max(diffractogram["2th"])-background_shoulders:
            peak_centers.append(new_peak)

    for peak_center in peak_centers:
        #print(peak_center)
        start_values_list.append(       [32.29992395, peak_center, 0.06344242,1])#[0.5,  peak_center,  0.1,1])#[0.1,  peak_center,  0.2434226])#,   1.0])
        peak_range_list.append(         [peak_center-peak_shoulders,peak_center+peak_shoulders])
        background_range_list.append(   [peak_center-background_shoulders,peak_center+background_shoulders])
    #print(start_values_list)
    fwhm_error_list=[]
    fwhm_parameter_list=[] 
    pos_error_list=[]
    pos_parameter_list=[]
    int_parameter_list=[]
    #using a fit to get values for fwhm in addition to intensity, peak_position and ratio of gaussian/lorentzian peak shape'''
    for j in range(len(peak_centers)):
        set_options= {
            'background_region': background_range_list[j],
            'peak_interval': peak_range_list[j],
            'background_excluded_region':False,
            'gaussian':False,
            #'voigt':False,
            'pseudovoigt': True,
            #Forcing the fit to be 100% gaussian by limiting the ratio between 0.9999 and 1
            'lower_bounds_PV': [0,peak_range_list[j][0],0,0.9999],#,0.99],#[0,peak_center-0.2,0,0,0],
            'upper_bounds_PV': [20000,peak_range_list[j][1],0.5,1]#,1.01],#[200,start_values_list[j][1]+0.2,0.5,1]
            #'plot_all_background_fits': file_to_plot
        }
        #Testing function by printing some of the peak fits
        #print(options)
        if j == 0:
            set_options['background_poly_degree']= 2
        else:
            set_options['background_poly_degree']= 1
        if j % options['plot_every_nth'] == 0:
            set_options['plot_fit']=True
            set_options['plot_all_background_fits'] = True
            
        else:
            set_options['plot_fit']=False
            set_options['plot_all_background_fits']: False

        output = background_subtracted_peak(data=data,options=set_options)
        df_peak=output[1]
        
        parameters, errors= find_fwhm_of_peak(x=df_peak["2th"],y=df_peak["I_corr"],start_values=start_values_list[j],options=set_options) 
        #Visualizing peak_center and fwhm
        #[I,x0,PV_fwhm,ratio]=parameters
        #if set_options['plot_fit']:
        #    plt.axvline(x0) #center position
        #    x=[x0-PV_fwhm/2,x0+PV_fwhm/2]
            #FIXME generalize for other functions than PV:
        #    sigma=PV_fwhm/(2*np.sqrt(2*np.log(2)))
        #    a_G= 1/(sigma*np.sqrt(2*np.pi))
    
            #b_G= 4*np.log(2)/PV_fwhm**2
            #LORENTZIAN_PART= 1/np.pi * (PV_fwhm/2)/((x_fit-x0)**2+(PV_fwhm/2)**2) #DISREGARDING THIS DUE TO PURELY GAUSSIAN
            #GAUSSIAN_PART= a_G*np.exp(-b_G*(x_fit-x0)**2)

        #y_fit = (I * (ratio * GAUSSIAN_PART+(1-ratio)*LORENTZIAN_PART))
        #    y=[(I*a_G)/2,(I*a_G)/2]
        #    plt.plot(x,y)

        
        int_parameter_list.append(parameters[0])

        pos_error_list.append(errors[1])
        pos_parameter_list.append(parameters[1])
        
        fwhm_error_list.append(errors[2])
        fwhm_parameter_list.append(parameters[2])
        #print(parameters)
    return int_parameter_list, pos_parameter_list, fwhm_parameter_list, fwhm_error_list

#def from_beamtime_to_wavelength(name_of_beamtime):
#    name_of_beamtime = name_of_beamtime[0].lower().replace('-', '').replace('_', '')
#    if name_of_beamtime == 'bm01021231':
#        print("beam time matches")
#        return 0.6390512
#    # Add more conditions for other beamtime names and wavelengths
#    else:
#        print("beam time does not match")
#        return None # Or raise an exception if the input name is invalid

def from_beamtime_to_wavelength(beamtime_name):
    # Convert the input to a list if it's not already
    if isinstance(beamtime_name, str):
        beamtime_name = [beamtime_name]
    # Clean up the beamtime names in the list
    beamtime_name = [name.lower().replace('-', '').replace('_', '') for name in beamtime_name]
    # Check if any of the beamtime names match and return the corresponding wavelength
    for name in beamtime_name:
        if name == 'bm01021231':
            return 0.6390512
        if name == 'bm01021239':
            return 0.6223090446548326
        if name == 'bm01021257':
            return 0.6892469162829824
        # Add more conditions for other beamtime names and wavelengths
    # If none of the beamtime names match, return None
    return None


def read_peak_width_from_refinement(filename):
    ad = None
    bd = None
    cd = None
    if isinstance(filename, list):
        filename=filename[0]
    with open(filename) as file:
        for line in file:
            if 'prm !ad' in line:
                ad = float(line.split('=')[1].strip('; \n'))
            elif 'prm !bd' in line:
                bd = float(line.split('=')[1].strip('; \n'))
            elif 'prm !cd' in line:
                cd = float(line.split('=')[1].strip('; \n'))
    return ad, bd, cd

import re

def get_bm_folder(path):
    ####################################### 
    # This is a function made by ChatGTP, that basically finds out what 
    # beamtime a dataset is from, if there exist a folder in the path 
    # starting with "bm" and contains 8 numbers (does not matteri if separated 
    # by dash/underline etc)
    ############################################
    # Split the path into its components
    components = os.path.normpath(path).split(os.sep)
    
    # Search for the first folder that matches the specified pattern
    pattern = r'^(bm|BM)\d{8}$'
    for component in components:
        match = re.match(pattern, component)
        if match:
            return match.group(0)
    
    # Return None if no matching folder was found
    return None

import os
import re

def get_bm_folder_v2(path):
    # Convert input to a list if it is not already a list
    if not isinstance(path, list):
        path = [path]

    bm_folders = []
    # Search for matching folders in each path in the input list
    for p in path:
        # Split the path into its components
        components = os.path.normpath(p).split(os.sep)
        
        # Search for the first folder that matches the specified pattern
        pattern = r'^(bm|BM)\d{8}$'
        for component in components:
            match = re.match(pattern, component)
            if match:
                bm_folders.append(match.group(0))
                break
            
    # Return the list of matching folders if the input was a list of paths
    if len(bm_folders) > 1:
        return bm_folders
    # Otherwise, return the name of the bm-folder as output
    elif len(bm_folders) == 1:
        return bm_folders[0]
    # If no matching folder was found, return None
    else:
        return None

def instrumental_peak_shape_plot(data,options):
    default_options = {
        'plot_instrumental_broadening': False,
        'refinement_result_path': None
    }
    
    options = aux.update_options(options=options, default_options=default_options)

    diffractogram, wavelength = xrd.io.read_xy(data=data,options=options)  

    ###### These functions can be put elsewhere, does not serve its purpose here? #####
    #beamtime = analyze.get_bm_folder_v2(data['path'])
    #wavelength = analyze.from_beamtime_to_wavelength(beamtime)
   #####################################################
    

    ad, bd, cd = read_peak_width_from_refinement(options['refinement_result_path'])
    min_x=min(diffractogram["2th"])
    max_x=max(diffractogram["2th"])

    x_fit=np.linspace(min_x,max_x,100)
    y_dmitry=_DC1(x_fit, ad, bd, cd)

    if options['plot_instrumental_broadening']:
        plt.plot(x_fit,y_dmitry)
        plt.title("LaB$_6$ as a function of angle")
        plt.xlabel("2th (deg)")
        plt.ylabel("fwhm")
    return ad,bd,cd

    #ChatGTP made this file (I though had to fix somethings myself) to find a file (with a certain string in the file name) in a folder and provide the full path of that file
def find_file(folder_path, detector_positions, file_string=None, file_exts=[".dat"]):
    matching_files = []
    for filename in os.listdir(folder_path):
        match_count = 0
        for detector_position in detector_positions:
            if detector_position in filename:
                match_count += 1
                break
        if file_string and file_string in filename:
            match_count += 1
        _, extension = os.path.splitext(filename)
        if extension in file_exts:
            match_count += 1
        if match_count == 3:
            matching_files.append(os.path.join(folder_path, filename))
    if len(matching_files) == 0:
        print(f"No files include the specified string(s) or extension(s) in folder {folder_path}")
        return None
    elif len(matching_files) == 1:
        return matching_files[0]
    else:
        print(f"Several files include the specified string(s) or extension(s) in folder {folder_path}")
        for file in matching_files:
            print("The LaB$_6$ xy-file used to define x min and max is "+str(file))
        return None











###Testing to make a function that only takes beam time as an argument:
def instrumental_peak_shape(beamtime,detector_position,options):
    default_options = {
        'plot_instrumental_broadening': False
        }

    options = aux.update_options(options=options, default_options=default_options)

    PATH=r"C:/Users/halvorhv/OneDriveUiO/0_Analysis_essentials"
    PATH_BEAMTIME=os.path.join(PATH,beamtime,"lab6")
    PATH_BEAMTIME_RESULTS=os.path.join(PATH_BEAMTIME,"results")
    #print(PATH_BEAMTIME)
    #print(PATH_BEAMTIME,["_000_00_","_100_80_"],"lab6","xy")
    if detector_position == "pos1":
        detector_positions = ["_000_00_","_100_80_","_000_080_"] #representing 1231 and 1257, more to come
    if detector_position == "pos3":
        detector_positions = ["_400_100_","_400_180_"] #representing 1257, more to come

    xy_path = find_file(PATH_BEAMTIME,detector_positions,"lab6",[".xy",".xye"])
    peak_width_path=find_file(PATH_BEAMTIME_RESULTS,detector_positions,"peak_width",[".dat"])

    if peak_width_path is None:
        print("Cannot find the Lab6 result-file, check detector position")
    
    wavelength= from_beamtime_to_wavelength(beamtime)
    

    data={
        'path': [xy_path],
        'wavelength': wavelength
    }

    diffractogram, wavelength = xrd.io.read_xy(data=data)#,options=options)  

    ###### These functions can be put elsewhere, does not serve its purpose here? #####
    #beamtime = analyze.get_bm_folder_v2(data['path'])
    #wavelength = analyze.from_beamtime_to_wavelength(beamtime)
   #####################################################
    
    ad, bd, cd = read_peak_width_from_refinement(peak_width_path)
    min_x=min(diffractogram["2th"])
    max_x=max(diffractogram["2th"])

    x_fit=np.linspace(min_x,max_x,100)
    y_dmitry=_DC1(x_fit, ad, bd, cd)

    if options['plot_instrumental_broadening']:
        plt.plot(x_fit,y_dmitry)
        plt.title("LaB$_6$ as a function of angle")
        plt.xlabel("2th (deg)")
        plt.ylabel("fwhm")
    return ad,bd,cd

def find_cbf_files(path):
    cbf_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".cbf"):
                cbf_files.append(os.path.join(root, file))
    return cbf_files

def from_lpa_Pt_to_temp(lpa_Pt):
    #Update a_0-value
    #FIXME add option of adding lpa_Pt_err
       #FIXME make new version with refining zero etc, and removing limits of scale_dis
#First finding the relation between relative expansion and temperature from Kirkby 1991
    T=[293.15,300,350,400,500,600,700,800,900,1000,1100,1200,1300]
    relativeexpansion=[0,61,512,971,1910,2871,3858,4870,5909,6976,8076,9213,10390]
    T1=285
    T2=480
    T3=1250
    T4=1800
    #z = np.polyfit(T, relativeexpansion, 3)
    x=np.linspace(T1,T4,num=T4-T1)
    #fit=z[0]*x*x*x+z[1]*x*x+z[2]*x+z[3]
    #print(fit)

    #ax=plt.plot(x, fit, '--', linewidth=2, markersize=12)

    #285-480 is recomended from Kirkby to be a separate fit
    Tlow=[293.15,300,350,400]; rel_exp_low=[0,61,512,971]
    Tmed=[500,600,700,800,900,1000,1100,1200]; rel_exp_med=[1910,2871,3858,4870,5909,6976,8076,9213]
    Thigh=[1300,1400,1500,1600,1700,1800] ;rel_exp_high=[10390,11604,12856,14151,15497,16908]

    zlow = np.polyfit(Tlow, rel_exp_low, 3)
    xlow=np.linspace(T1,T2,num=T2-T1)
    #FITTING THE LOW TEMPERATURE REGION
    fit_low=zlow[0]*xlow*xlow*xlow+zlow[1]*xlow*xlow+zlow[2]*xlow+zlow[3]

    #bx=plt.plot(Tlow, rel_exp_low, 'go', linewidth=2, markersize=5)
    #bx=plt.plot(xlow, fit_low, '--', linewidth=2, markersize=12)
    zmed = np.polyfit(Tmed, rel_exp_med, 3)
    xmed=np.linspace(T2,T3,num=T3-T2)
    #FITTING HIGH TEMPERATURE REGION
    fit_med=zmed[0]*xmed*xmed*xmed+zmed[1]*xmed*xmed+zmed[2]*xmed+zmed[3]
    #print(fit_med)
    zhigh = np.polyfit(Thigh, rel_exp_high, 3)
    xhigh=np.linspace(T3,T4,num=T4-T3)
    #FITTING HIGH TEMPERATURE REGION
    fit_high=zhigh[0]*xhigh*xhigh*xhigh+zhigh[1]*xhigh*xhigh+zhigh[2]*xhigh+zhigh[3]

    fit_new1=np.append(fit_low,fit_med)
    fit_new=np.append(fit_new1,fit_high)

    #PLOTTING THE DATA POINTS
    T=np.array(T)# - 293
    x=np.array(x)
    T_C=T-273.15
    x_C=x-273.15

    #T_C = [T - 293 for x in T]

    #IF plotting:
    #ax=plt.plot(T_C, relativeexpansion, 'go', linewidth=2, markersize=5)

    #bx=plt.plot(x_C, fit_new, '--', linewidth=2, markersize=12)

    ## ==== value below gave good fit with blower temperatures ...
    #a_0 = float(3.926000666666667) #from 2023_temperature_calibration.ipynb
    #a_0_error = float(4.4743714642394187e-05) #from 2023_temperature_calibration.ipynb
  
    a_0_Mn15_8Q_Pt_02=3.922855 #fra v9
    a_0_Mn15_8Q_Pt_01=3.923074 #fra v9
    a_0_Mn15_8Q_O2_slowcool_Pt=3.922450 #fra v9
    a_0 = float(np.average([a_0_Mn15_8Q_Pt_02,a_0_Mn15_8Q_Pt_01,a_0_Mn15_8Q_O2_slowcool_Pt]))#3.929536666666667) #too high????

    
    relative_expansion = (lpa_Pt - a_0)*1000000/a_0

    # using enumerate() + next() to find index of first element in fit_new (temperature calibration) just greater than the calculated halfmax 
    realtemp_index = next(x for x, val in enumerate(fit_new) 
                        if val > relative_expansion)
    
    #print(realtemp_index)
    temp = x_C[realtemp_index]
    return temp

def make_column_in_df_with_calib_temp(df,options):
# dataframe needs a column with the name "lpa_Pt" for this to work.
    default_options = {
        'plot': False
        }

    options = aux.update_options(options=options, default_options=default_options)

    calib_temp=[]
    #Plan: Find which row in newnames-column of temp_new.txt contains the "filenumber" from each row of df_new. Take the blower-value from this row and add to a new column in df_new: "blowertemp".
    for index, row in df.iterrows(): #iterating through each row  
        lpa_Pt = row['lpa_Pt']
        temp = from_lpa_Pt_to_temp(lpa_Pt)
        calib_temp.append(temp)
        
    df["Calib_temp"]=calib_temp
    if options["plot"]:
        df.plot(y=['Blower','Calib_temp'], use_index=True)
    return df

def line_prepender(filename, line):
        #function that adds a line in the beginning of a file (adapted for big.inps)
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

import os

def line_insert_after(filename, line,string_to_locate):
    #function to add a line of code after a certain string in the file. Ideal for sequential refinement, adding filenames to the .inp-file.
    temp_filename = filename + '.tmp'
    with open(filename, 'r') as f:
        with open(temp_filename, 'w') as temp_file:
            for current_line in f:
                temp_file.write(current_line)
                if string_to_locate in current_line:
                    temp_file.write(line.rstrip('\r\n') + '\n')

    os.remove(filename)
    os.rename(temp_filename, filename)

def find_closest_filename(df, target_time):
# Calculate the absolute time difference between 'Rel_time' and the target_time in in situ sequential analysis
    df['time_difference'] = abs(df['Rel_time'] - target_time)

    # Sort the DataFrame based on the 'time_difference' column in ascending order
    df_sorted = df.sort_values('time_difference')

    # Get the filename from the row with the smallest 'time_difference'
    closest_filename = df_sorted.iloc[0]['filename']

    # Drop the temporary 'time_difference' column
    df.drop('time_difference', axis=1, inplace=True)

    return closest_filename

def set_value_as_macro(s, parameter_name,macro_name):
    string = s.rsplit("local "+parameter_name)[-1]

    value_to_replace = string.split("min")[0]
    print(value_to_replace)
    s = s.replace(value_to_replace,"= "+macro_name+" ;:	 0` ")
    s = s.replace("local "+parameter_name,"local !"+parameter_name)
    return s