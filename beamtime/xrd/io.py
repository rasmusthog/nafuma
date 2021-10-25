import fabio, pyFAI
import pandas as pd
import numpy as np
import os

import zipfile
import io
import sys
import glob
import xml.etree.ElementTree as ET


def get_image_array(path):

    image = fabio.open(path)
    image_array = image.data

    return image_array


def integrate_1d(calibrant, bins, path=None, image=None, options=None):
    ''' Integrates an image file to a 1D diffractogram. 

    Input:
    calibrant: path to .poni-file
    bins: Number of bins to divide image into
    path (optional): path to image file - either this or image must be specified. If both is passed, image is prioritsed
    image (optional): image array (Numpy) as extracted from get_image_array
    options (optional): dictionary of options

    Output:
    df: DataFrame contianing 1D diffractogram if option 'return' is True
    ''' 

    required_options = ['unit', 'extension', 'filename', 'save_folder', 'overwrite', 'return']

    default_options = {
        'unit': '2th_deg', 
        'extension': '_integrated.dat',
        'filename': None,
        'save_folder': '.',
        'overwrite': False,
        'return': False}

    if not options:
        options = default_options
    
    else:
        for option in required_options:
            if option not in options.keys():
                options[option] = default_options[option]
    
    
    if not image:
        image = get_image_array(path)
    
    ai = pyFAI.load(calibrant)


    if not options['filename']:
        if path:
            filename = os.path.join(options['save_folder'], os.path.split(path)[-1].split('.')[0] + options['extension'])
        else:
            filename = os.path.join(options['save_folder'], 'integrated.dat')


    if not options['overwrite']:
        trunk = os.path.join(options['save_folder'], filename.split('\\')[-1].split('.')[0])
        extension = filename.split('.')[-1]
        counter = 0

        while os.path.isfile(filename):
            counter_string = str(counter)
            filename = trunk + '_' + counter_string.zfill(4) + '.' + extension
            counter += 1


    if not os.path.isdir(options['save_folder']):
        os.makedirs(options['save_folder'])


    res = ai.integrate1d(image, bins, unit=options['unit'], filename=filename)

    if options['return']:
        return open_1d_data(filename)
    

def open_1d_data(path, options=None):

    with open(path, 'r') as f:
        position = 0
    
        current_line = f.readline()
        
        while current_line[0] == '#':
            position = f.tell()
            current_line = f.readline()
            
        f.seek(position)

        df = pd.read_csv(f, header=None, delim_whitespace=True)

    df.columns = ['2th', 'I']


    return df
    


def average_images(images):
    ''' Takes a list of path to image files, reads them and averages them before returning the average image'''


    image_arrays = []

    for image in images:
        image_array = xrd.io.get_image_array(os.path.join(root, image))
        image_arrays.append(image_array)
    
    
    image_arrays = np.array(image_arrays)

    image_average = image_arrays.mean(axis=0)


    return image_average


def subtract_dark(image, dark):

    return image - dark



def view_integrator(calibrant):
    ''' Prints out information about the azimuthal integrator
    
    Input:
    calibrant: Path to the azimuthal integrator file (.PONI)
    
    Output:
    None'''

    ai = pyFAI.load(calibrant)

    print("pyFAI version:", pyFAI.version)
    print("\nIntegrator: \n", ai)


def brml_reader(file_name):
    """
    From: https://github.com/aboulle/DxTools/blob/master/data_reader.py
    
    Extracts brml into a temporary unzip file, and parses all xml files.
    For every intensity value, saves all relevant motor coordinates and sensor values.
    All values are save in temporary SPEC-style tmp file.
    """
#*****************************************************************************************************
# Unzip xml files
#*****************************************************************************************************
    extract_path = os.path.join(os.getcwd(),"unzip")
    if sys.platform == "win32":
        print("Detected platform: Windows")
        os.system("RMDIR "+ extract_path +" /s /q")
    elif sys.platform == "darwin":
        print("Detected platform: MacOS")
        os.system("rm -rf "+ extract_path)
    elif sys.platform == "linux" or sys.platform == "linux2":
        print("Detected platform: Linux")
        os.system("rm -rf "+ extract_path)

    #Extract all RawData*.xml files and InstructionContainer the brml to temporary unzip file
    with zipfile.ZipFile(file_name,"r") as brml:
        for info in brml.infolist():
            if ("RawData" in info.filename) or ("InstructionContainer" in info.filename):
            #if ("RawData" in info.filename):
                brml.extract(info.filename, extract_path)
#*****************************************************************************************************
# For time counting, the number of days is initialized to 0.
# Compatibility fixes with D8 advance and older Discover: offsets, chi and tx, ty are initialized to 0
#*****************************************************************************************************
    # Initialize the number of days to 0
    n_day = 0.
    # Initialize all offsets to 0
    off_tth = off_om = off_phi = off_chi = off_tx = off_ty = 0
    # Set Chi, tx and ty to 0 for D8 advance (July 2017 Julia Stroh)
    chi = tx = ty = "0"
#*****************************************************************************************************
#Modification June 2017 (Duc Dinh)
#In some RawData.xml files, wavelength and static motors are missing.
#Find wl and static motors in MeasurementContainer.xml
#*****************************************************************************************************
    data_path = os.path.join(extract_path, "*0","InstructionContainer.xml")
    for file in sorted(glob.glob(data_path)):
        tree = ET.parse(file)
        root = tree.getroot()
        for chain in root.findall("./ComparisonMethod/HrxrdAlignmentData"):
            wl = chain.find("WaveLength").attrib["Value"]

        for chain in root.findall("./ComparisonMethod/HrxrdAlignmentData/Data"):
            if chain.get("LogicName") == "TwoTheta":
                tth = chain.find("TheoreticalPosition").attrib["Value"]
                off_tth = chain.find("PositionOffset").attrib["Value"]
                tth = str(float(tth)-float(off_tth))

            if chain.get("LogicName") == "Theta":
                om = chain.find("TheoreticalPosition").attrib["Value"]
                off_om = chain.find("PositionOffset").attrib["Value"]
                om = str(float(om)-float(off_om))

            if chain.get("LogicName") == "Chi":
                chi = chain.find("TheoreticalPosition").attrib["Value"]
                off_chi = chain.find("PositionOffset").attrib["Value"]
                chi = str(float(chi)-float(off_chi))

            if chain.get("LogicName") == "Phi":
                phi = chain.find("TheoreticalPosition").attrib["Value"]
                off_phi = chain.find("PositionOffset").attrib["Value"]
                phi = str(float(phi)-float(off_phi))

            if chain.get("LogicName") == "X":
                tx = chain.find("TheoreticalPosition").attrib["Value"]
                off_tx = chain.find("PositionOffset").attrib["Value"]
                tx = str(float(tx)-float(off_tx))

            if chain.get("LogicName") == "Y":
                ty = chain.find("TheoreticalPosition").attrib["Value"]
                off_ty = chain.find("PositionOffset").attrib["Value"]
                ty = str(float(ty)-float(off_ty))
        os.remove(file)

    #Create ouput file
    outfile = open("tmp", "w", encoding='utf8') # Create output data file
    outfile.write("#temperature   khi   phi   x   y   theta   offset   2theta   scanning motor   intensity   time\n")
#*****************************************************************************************************
# Finds scan type, wl, scanning motors and fixed motors values in RawData*.xml
#*****************************************************************************************************
    data_path = os.path.join(extract_path, "*0","*.xml") #reading files in Experiment0 folder
    for file in sorted(glob.glob(data_path), key=file_nb):
        new_file = 0
        check_temperature = 0
        check_1Dmode = 0
        #parsing XML file
        tree = ET.parse(file) 
        root = tree.getroot()
        #obtain scan type
        for chain in (root.findall("./DataRoutes/DataRoute/ScanInformation") or root.findall("./ScanInformation")):
            scan_type = chain.get("VisibleName")
            if ("PSD" in scan_type) or ("Psd" in scan_type):
                scan_type = "PSDFIXED"
            if ("Coupled" in scan_type) or ("coupled" in scan_type) or ("2Theta-Omega" in scan_type):
                scan_type = "COUPLED"
            if ("Rocking" in scan_type) or ("rocking" in scan_type):
                scan_type = "THETA"

        # Check if temperature is recorded
        for chain in (root.findall("./DataRoutes/DataRoute/DataViews/RawDataView/Recording") or root.findall("./DataViews/RawDataView/Recording")):
            if "Temperature" in chain.get("LogicName"):
                check_temperature = 1

        #Find wl in RawData.xml
        for chain in root.findall("./FixedInformation/Instrument/PrimaryTracks/TrackInfoData/MountedOptics/InfoData/Tube/WaveLengthAlpha1"):
            wl = chain.get("Value")

        # Find the fast-scanning axis
        for chain in (root.findall("./DataRoutes/DataRoute/ScanInformation/ScanAxes/ScanAxisInfo") or root.findall("./ScanInformation/ScanAxes/ScanAxisInfo")):
            if new_file == 0:
                if chain.get("AxisName") == "TwoTheta": #Added offset correction / June 2017. Only relevant if offset in InstructionContainer. 0 otherwise.
                    off_scan = float(off_tth)
                elif chain.get("AxisName") == "Theta":
                    off_scan = float(off_om)
                else:
                    off_scan = 0
                step = chain.find("Increment").text
                start = chain.find("Start").text
                stop = chain.find("Stop").text
                ref = chain.find("Reference").text
                start = str(float(ref)+float(start)-off_scan)  #Added offset correction / June 2017.
                #start = str(float(ref)+float(start))
                new_file += 1

        # Find scanning motors
        for chain in (root.findall("./DataRoutes/DataRoute/ScanInformation/ScanAxes/ScanAxisInfo") or root.findall("./ScanInformation/ScanAxes/ScanAxisInfo")):
            if chain.get("AxisName") == "TwoTheta":
                tth = chain.find("Start").text
                ref = chain.find("Reference").text
                tth = str(float(ref)+float(tth)-float(off_tth))  #Added offset correction / June 2017.
                #tth = str(float(ref)+float(tth))

            if chain.get("AxisName") == "Theta":
                om = chain.find("Start").text
                ref = chain.find("Reference").text
                om = str(float(ref)+float(om)-float(off_om))  #Added offset correction / June 2017.
                #om = str(float(ref)+float(om))

            if chain.get("AxisName") == "Chi":
                chi = chain.find("Start").text
                ref = chain.find("Reference").text
                chi = str(float(ref)+float(chi)-float(off_chi))  #Added offset correction / June 2017.
                #chi = str(float(ref)+float(chi))

            if chain.get("AxisName") == "Phi":
                phi = chain.find("Start").text
                ref = chain.find("Reference").text
                phi = str(float(ref)+float(phi)-float(off_phi))  #Added offset correction / June 2017.
                #phi = str(float(ref)+float(phi))

            if chain.get("AxisName") == "X":
                tx = chain.find("Start").text
                ref = chain.find("Reference").text
                tx = str(float(ref)+float(tx)-float(off_tx))  #Added offset correction / June 2017.
                #tx = str(float(ref)+float(tx))

            if chain.get("AxisName") == "Y":
                ty = chain.find("Start").text
                ref = chain.find("Reference").text
                ty = str(float(ref)+float(ty)-float(off_ty))  #Added offset correction / June 2017.
                #ty = str(float(ref)+float(ty))

        # Find static motors
        for chain in root.findall("./FixedInformation/Drives/InfoData"):
            if chain.get("LogicName") == "TwoTheta":
                tth = chain.find("Position").attrib["Value"]
                tth = str(float(tth)-float(off_tth))  #Added offset correction / June 2017.

            if chain.get("LogicName") == "Theta":
                om = chain.find("Position").attrib["Value"]
                om = str(float(om)-float(off_om))  #Added offset correction / June 2017.

            if chain.get("LogicName") == "Chi":
                chi = chain.find("Position").attrib["Value"]
                chi = str(float(chi)-float(off_chi))  #Added offset correction / June 2017.

            if chain.get("LogicName") == "Phi":
                phi = chain.find("Position").attrib["Value"]
                phi = str(float(phi)-float(off_phi))  #Added offset correction / June 2017.

            if chain.get("LogicName") == "X":
                tx = chain.find("Position").attrib["Value"]
                tx = str(float(tx)-float(off_tx))  #Added offset correction / June 2017.

            if chain.get("LogicName") == "Y":
                ty = chain.find("Position").attrib["Value"]
                ty = str(float(ty)-float(off_ty))  #Added offset correction / June 2017.

        offset = str(float(om) - float(tth)/2.)

#*****************************************************************************************************
# This section computes scanning time, scanning angular range and scanning speed
# in order to convert 2th values to time values (July 2017, Julia Stroh)
#*****************************************************************************************************
        for chain in (root.findall("./TimeStampStarted")):
            d_start = ((chain.text).split("T")[0]).split("-")[2]
            t_start = ((chain.text).split("T")[1]).split("+")[0]
            h_start, min_start, sec_start = t_start.split(":")
            t_start = float(h_start)*3600 + float(min_start)*60 + float(sec_start)

            if file_nb(file)==0:
                abs_start = t_start
        for chain in (root.findall("./TimeStampFinished")):
            d_stop = ((chain.text).split("T")[0]).split("-")[2]
            t_stop = ((chain.text).split("T")[1]).split("+")[0]
            h_stop, min_stop, sec_stop = t_stop.split(":")
            t_stop = float(h_stop)*3600 + float(min_stop)*60 + float(sec_stop)

        # Check if detector is in 1D mode
        for chain in (root.findall("./DataRoutes/DataRoute/ScanInformation") or root.findall("./ScanInformation")):
            if chain.find("TimePerStep").text != chain.find("TimePerStepEffective").text:
                check_1Dmode = 1

        # Check if day changed between start and stop and correct accordingly
        if d_stop != d_start:
            t_stop += 24*3600.
        total_scan_time = t_stop - t_start

        #scanning range
        dth_scan = float(stop)-float(start)
        #psd range
        dth_psd = 0
        for chain in root.findall("./FixedInformation/Detectors/InfoData/AngularOpening"):
            dth_psd = chain.get("Value")
        total_dth = float(dth_psd)*check_1Dmode+float(dth_scan)

        scan_speed = total_dth / total_scan_time
#*****************************************************************************************************
# Finds intensity values. If temperature is recorded, also fin temperature values.
# The intensity data is formatted differently in PSDfixed mode and when temperature is recorded
#*****************************************************************************************************
        if "PSDFIXED" in scan_type:
            if check_temperature == 0:
                for chain in (root.findall("./DataRoutes/DataRoute") or root.findall("./")):
                    intensity = (chain.find("Datum").text).split(',')

                for chain in (root.findall("./DataRoutes/DataRoute/DataViews/RawDataView/Recording") or root.findall("./DataViews/RawDataView/Recording")):
                    if chain.get("LogicName") == "Counter1D":
                        n_channels = int(chain.find("Size/X").text)

                line_count = 0
                int_shift = len(intensity) - n_channels
                for i in range(n_channels): #the intensity values are shifted to the right by int_shift
                    if i == 0:
                        scanning = float(start)
                    else:
                        scanning += float(step)
                    line_count += 1
                    t_2th = (t_start+n_day*24*3600 - abs_start)+((float(dth_psd)*check_1Dmode + scanning - float(start)) / scan_speed)
                    outfile.write("25" + " " + (chi) + " " + (phi)
                                  + " " + (tx) + " " + (ty) + " " + (om)
                                  + " " + (offset) + " " + (tth) + " " + str(scanning)
                                  + " "  + intensity[i+int_shift] +" " + str(t_2th) +'\n')
            else:
                return implementation_warning, 0, 0

        #if "COUPLED" in scan_type:
        # to do check in brml that all scans (except psd fixed) share the same data structure (wrt temperature)
        else:
            if check_temperature == 0:
                line_count = 0
                for chain in (root.findall("./DataRoutes/DataRoute/Datum") or root.findall("./Datum")):
                    if line_count == 0:
                        scanning = float(start)
                    else:
                        scanning += float(step)
                    line_count += 1
                    intensity = (chain.text).split(',')[-1]
                    #compute time corresponding to scanning angle (July 2017)
                    t_2th = (t_start+n_day*24*3600 - abs_start)+((float(dth_psd)*check_1Dmode + scanning - float(start)) / scan_speed)
                    outfile.write("25" + " " + (chi) + " " + (phi)
                                  + " " + (tx) + " " + (ty) + " " + (om)
                                  + " " + (offset) + " " + (tth) + " " + str(round(scanning, 4))
                                  + " " + intensity + " " + str(t_2th) +'\n')
            else:
                line_count = 0
                for chain in (root.findall("./DataRoutes/DataRoute/Datum") or root.findall("./Datum")):
                    if line_count == 0:
                        scanning = float(start)
                    else:
                        scanning += float(step)
                    line_count += 1
                    t_2th = (t_start+n_day*24*3600 - abs_start)+((float(dth_psd)*check_1Dmode + scanning - float(start)) / scan_speed)
                    intensity = (chain.text).split(',')[-2]
                    temperature = (chain.text).split(',')[-1]
                    outfile.write(temperature + " " + (chi) + " " + (phi)
                                  + " " + (tx) + " " + (ty) + " " + (om)
                                  + " " + (offset) + " " + (tth) + " " + str(round(scanning, 4))
                                  + " "  + intensity + " " + str(t_2th) +'\n')

        if d_stop != d_start:
            n_day+=1
            
    outfile.close()
    if sys.platform == "win32":
        os.system("RMDIR "+ extract_path +" /s /q")
    elif sys.platform == "darwin":
        os.system("rm -rf "+ extract_path)
    elif sys.platform == "linux" or sys.platform == "linux2":
        os.system("rm -rf "+ extract_path)
    return scan_type, line_count, wl