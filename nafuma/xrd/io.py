from sympy import re
import fabio, pyFAI
import pandas as pd
import numpy as np
import os
import shutil
import sys
import datetime

import zipfile
import xml.etree.ElementTree as ET


import nafuma.auxillary as aux


def get_image_array(path):

    beamline_extension = ['.edf', '.cbf', '.mar3450']

    if path.endswith(tuple(beamline_extension)):
        image = fabio.open(path)
        image_array = image.data

    elif path.endswith('.dat'):
        image_array = np.loadtxt(path, skiprows=1, delimiter=';')

    return image_array


def get_image_headers(path):

    image = fabio.open(path)

    return image.header



def integrate_scans(data: dict, options={}):

    default_options = {
        'extension': '.dat',
        'save': True,
        'integration_save_folder': './integrated/',
        'filename_base': 'integrated',
    }

    options = aux.update_options(options=options, required_options=default_options.keys(), default_options=default_options)


    if not isinstance(data['path'], list):
        imgs = aux.get_filenames(data['path'], ext=options['extension'])


    diffractograms, wavelengths = [], []
    for i, img in enumerate(imgs):
        data['image'] = get_image_array(img)

        options['integration_save_filename'] = options['filename_base'] + '_' + f'{i}'.zfill(4) + '.xy'

        diff, wl = integrate_1d(data=data, options=options)

        diffractograms.append(diff)
        wavelengths.append(wl)
    
    return diffractograms, wavelengths


def integrate_1d(data, options={}, index=0):
    ''' Integrates an image file to a 1D diffractogram. 

    Required content of data:
    calibrant (str): path to .poni-file
    nbins (int): Number of bins to divide image into
    path (str) (optional, dependent on image): path to image file - either this or image must be specified. If both is passed, image is prioritsed
    image (NumPy 2D Array) (optional, dependent on path): image array as extracted from get_image_array

    Output:
    df: DataFrame contianing 1D diffractogram if option 'return' is True
    ''' 

    required_options = ['unit', 'npt', 'save', 'integration_save_filename', 'save_extension', 'integration_save_folder', 'overwrite', 'extract_folder', 'error_model']

    default_options = {
        'unit': '2th_deg', 
        'npt': 5000,
        'extract_folder': 'tmp',
        'error_model': None,
        'save': False,
        'integration_save_filename': None,
        'save_extension': '_integrated.xy',
        'integration_save_folder': '.',
        'overwrite': False}

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]
    

    # Get image array from filename if not passed
    if 'image' not in data.keys() or not isinstance(data['image'], np.ndarray):
        data['image'] = get_image_array(data['path'][index])


    # Load mask
    if 'mask' in data.keys():
        mask = get_image_array(data['mask'])
    else:
        mask = None
    

    # Instanciate the azimuthal integrator from pyFAI from the calibrant (.poni-file)
    ai = pyFAI.load(data['calibrant'])

    # Determine filename
    filename = make_filename(options=options, path=data['path'][index])
    
    # Make save_folder if this does not exist already
    if not os.path.isdir(options['extract_folder']):
        os.makedirs(options['extract_folder'])

    if not os.path.isdir(options['integration_save_folder']):
        os.makedirs(options['integration_save_folder'])

    


    res = ai.integrate1d(data['image'], npt=options['npt'], mask=mask, error_model=options['error_model'], unit=options['unit'], filename=filename)

    data['path'][index] = filename
    diffractogram, _ = read_xy(data=data, options=options, index=index)
    wavelength = find_wavelength_from_poni(path=data['calibrant'])

    if not options['save']:
        os.remove(filename)
        shutil.rmtree(f'tmp')
    
    return diffractogram, wavelength
    


def make_filename(options, path=None):

    # Define save location for integrated diffractogram data
    if not options['save']:
            filename = os.path.join(options['extract_folder'], 'tmp_diff.dat')

    elif options['save']:

        # Case 1: No filename is given. 
        if not options['integration_save_filename']:
            # If a path is given instead of an image array, the path is taken as the trunk of the savename
            if path:
                # Make filename by joining the save_folder, the filename (with extension deleted) and adding the save_extension
                filename = os.path.join(options['integration_save_folder'], os.path.split(path)[-1].split('.')[0] + options['save_extension'])
            else:
                # Make filename just "integrated.dat" in the save_folder
                filename = os.path.join(options['integration_save_folder'], 'integrated.xy')


        else:
            filename = os.path.join(options['integration_save_folder'], options['integration_save_filename'])

        if not options['overwrite']:
            trunk = filename.split('.')[0]
            extension = filename.split('.')[-1]
            counter = 0

            while os.path.isfile(filename):
                
                # Rename first file to match naming scheme if already exists
                if counter == 0:
                    os.rename(filename, trunk + '_' + str(counter).zfill(4) + '.' + extension)
                
                # Increment counter and make new filename
                counter += 1
                counter_string = str(counter)
                filename = trunk + '_' + counter_string.zfill(4) + '.' + extension


    return filename
    


def generate_image_list(path, options=None):
    ''' Generates a list of paths to pass to the average_images() function'''

    required_options = ['scans_per_image']
    default_options = {
        'scans_per_image': 5
    }


def process_2d_scans(data: dict, options={}):

    default_options = {
        'scans': 15, # number of scans per image
        'img_filename': 'img_',
        'extension': '.edf',
        'darks': True, # whether there are darks
        'dark_filename': 'dark_',
        'save': False,
        'save_folder': './average/',
        'save_filename': 'avg_',
        'save_extension': '.dat'
    }

    options = aux.update_options(options=options, required_options=default_options.keys(), default_options=default_options)


    all_imgs = [os.path.join(data['path'], img) for img in os.listdir(data['path']) if img.endswith(options['extension']) and img.startswith(options['img_filename'])]

    if options['darks']:
        all_darks = [os.path.join(data['path'], img) for img in os.listdir(data['path']) if img.endswith(options['extension']) and img.startswith(options['dark_filename'])]


    scans = int(len(all_imgs) / options['scans'])


    assert scans - (len(all_imgs) / options['scans']) == 0


    imgs = []
    darks = []

    for i in range(scans):
        img = []
        dark = []

        for j in range(options['scans']):
            img.append(all_imgs.pop(0))

            if options['darks']:
                dark.append(all_darks.pop(0))

        imgs.append(img)
        
        if options['darks']:
            darks.append(dark)


    img_avgs = [] 
    headers = []

    for img, dark in zip(imgs,darks):
        img_avg = average_images(img)
        header = get_image_headers(img[0])  
        
        if options['darks']:
            dark_avg = average_images(dark)
            img_avg = subtract_dark(img_avg, dark_avg)

        img_avgs.append(img_avg)
        headers.append(header)


    if options['save']:
        if not os.path.isdir(options['save_folder']):
            os.makedirs(options['save_folder'])

        for i, img in enumerate(img_avgs):
            if options['save_extension'] == '.dat':
                with open(os.path.join(options['save_folder'], options['save_filename']+f'{i}'.zfill(4)+options['save_extension']), 'w') as f:
                    f.write(f'# Time: {headers[i]["time"]}\n')
                    np.savetxt(f, img, fmt='%.2f', delimiter=";")
                

    return img_avgs



        





def average_images(images):
    ''' Takes a list of path to image files, reads them and averages them before returning the average image'''

    image_arrays = []

    for image in images:
        image_array = get_image_array(image)
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




def read_brml(data, options={}, index=0):


    # FIXME: Can't read RECX1-data, apparently must be formatted differently from RECX2. Check the RawData-files and compare between the two files.


    required_options = ['extract_folder', 'save_folder']
    default_options = {
        'extract_folder': 'tmp',
        'save_folder': None
    }


    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    if not os.path.isdir(options['extract_folder']):
        os.mkdir(options['extract_folder'])


    # Extract the RawData0.xml file from the brml-file
    with zipfile.ZipFile(data['path'][index], 'r') as brml:
        for info in brml.infolist():
            if "RawData" in info.filename:
                brml.extract(info.filename, options['extract_folder'])



    # Parse the RawData0.xml file
    path = os.path.join(options['extract_folder'], 'Experiment0/RawData0.xml')

    tree = ET.parse(path)
    root = tree.getroot()

    shutil.rmtree(options['extract_folder'])

    diffractogram = []

    for chain in root.findall('./DataRoutes/DataRoute'):


        # Get the scan type to be able to handle different data formats
        scantype = chain.findall('ScanInformation')[0].get('VisibleName')

        # Check if the chain is the right one to extract the data from
        if chain.get('Description') == 'Originally measured data.':
            
            
            if scantype == 'TwoTheta':
                for scandata in chain.findall('Datum'):
                    scandata = scandata.text.split(',')
                    twotheta, intensity = float(scandata[2]), float(scandata[3])
                    
                    if twotheta > 0:
                        diffractogram.append({'2th': twotheta, 'I': intensity})

            elif scantype == 'Coupled TwoTheta/Theta':
                for scandata in chain.findall('Datum'):
                    scandata = scandata.text.split(',')
                    twotheta, intensity = float(scandata[2]), float(scandata[4])
                    
                    if twotheta > 0:
                        diffractogram.append({'2th': twotheta, 'I': intensity})

            elif scantype == 'Still (Eiger2R_500K (1D mode))':

                start = float(chain.findall('ScanInformation/ScaleAxes/ScaleAxisInfo/Start')[0].text)
                stop = float(chain.findall('ScanInformation/ScaleAxes/ScaleAxisInfo/Stop')[0].text)
                



                for scandata in chain.findall('Datum'):
                        scandata = scandata.text.split(',')
                        raw = [float(i) for i in scandata]

                        intensity = []
                        for r in raw:
                            if r > 601:
                                intensity.append(r)

                        intensity = np.array(intensity)

            
             

                twotheta = np.linspace(start, stop, len(intensity))

                diffractogram = {'2th': twotheta, 'I': intensity}



    #if 'wavelength' not in data.keys():
    # Find wavelength
    
    if not data['wavelength'][index]:
        for chain in root.findall('./FixedInformation/Instrument/PrimaryTracks/TrackInfoData/MountedOptics/InfoData/Tube/WaveLengthAlpha1'):
            wavelength = float(chain.attrib['Value'])
    else:
        wavelength = data['wavelength'][index]


    diffractogram = pd.DataFrame(diffractogram)
    



    if options['save_folder']:
        if not os.path.isdir(options['save_folder']):
            os.makedirs(options['save_folder'])

        diffractogram.to_csv(options['save_folder'])



    return diffractogram, wavelength
    

def read_htxrd(data, options={}, index=0):

    required_options = ['extract_folder', 'save_folder', 'save_filename', 'adjust_time']
    default_options = {
        'extract_folder': 'tmp',
        'save_folder': None,
        'save_filename': None,
        'adjust_time': True
    }

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

    if 'wavelength' not in data.keys():
        data['wavelength'] = [None for i in range(len(data['path']))]


    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    
    # Extract the RawData0.xml file from the brml-file
    with zipfile.ZipFile(data['path'][index], 'r') as brml:
        for info in brml.infolist():
            if "RawData" in info.filename:
                brml.extract(info.filename, options['extract_folder'])


    # Get all filenames
    files = os.listdir(os.path.join(options['extract_folder'], 'Experiment0'))

    # initalise empty list to store all DataFrames
    diffractograms = []
    wavelengths = []
    
    active_scan = False
    timestamps = []

    # Loop through all RawData-files and extract all data and temperatures
    for i, file in enumerate(files):

        # Create all filenames as strings
        filename = os.path.join('tmp/Experiment0/', f'RawData{i}.xml')

        # Parse the .xml-files
        tree = ET.parse(filename)
        root = tree.getroot()

        # initalise empty list to store data from this particular scan
        diffractogram = []

        for chain in root.findall('./DataRoutes/DataRoute'):
            
            scantypes = chain.findall('ScanInformation')

            for scantype in scantypes:
                if scantype.get('VisibleName') == 'Still (TCU1000N)':
                    continue

                else:
                    active_scan = True
                    if chain.get('RouteFlag') == 'Final':
                        for scandata in chain.findall('Datum'):
                            scandata = scandata.text.split(',')
                            twotheta, intensity, temperature = float(scandata[2]), float(scandata[3]), float(scandata[5])

                            diffractogram.append({'2th': twotheta, 'I': intensity, 'T': temperature})

                        diffractogram = pd.DataFrame(diffractogram)
                        diffractograms.append(diffractogram)


                        if not data['wavelength'][index]:
                            for chain in root.findall('./FixedInformation/Instrument/PrimaryTracks/TrackInfoData/MountedOptics/InfoData/Tube/WaveLengthAlpha1'):
                                wavelength = float(chain.attrib['Value'])
                        else:
                            wavelength = data['wavelength'][index]

                        wavelengths.append(wavelength)


        if active_scan:
            for chain in root.findall('./TimeStampStarted'):
                time_start = datetime.datetime.strptime(chain.text[:-7], "%Y-%m-%dT%H:%M:%S.%f")
            for chain in root.findall('./TimeStampFinished'):
                time_end = datetime.datetime.strptime(chain.text[:-7], "%Y-%m-%dT%H:%M:%S.%f")


            time_diff = time_end - time_start

            if options['adjust_time']:
                timestamps.append(time_start + time_diff/2)


    if options['save_folder']:
        for i, (diffractogram, wavelength, timestamp) in enumerate(zip(diffractograms, wavelengths, timestamps)):
            if not options['save_filename']:
                filename = os.path.basename(data['path'][index]).split('.')[0] + '_' + str(i).zfill(4) +'.xy'
            else:
                filename = options['save_filename'] + '_' + str(i).zfill(4) +'.xy'


            if not os.path.isdir(options['save_folder']):
                os.makedirs(options['save_folder'])

            save_htxrd_as_xy(diffractogram, wavelength, timestamp, filename, options['save_folder'])





    shutil.rmtree(options['extract_folder'])

    return diffractograms, wavelengths

def save_htxrd_as_xy(diffractogram, wavelength, timestamp, filename, save_path):

    headers = '\n'.join(
        [line for line in 
            [f'# Temperature {np.round(diffractogram["T"].mean())}',
            f'# Wavelength {wavelength}',
            f'# Time {timestamp}',
            '# 2th \t I'
            ]
        ]
    )

    diffractogram = diffractogram.drop('T', axis=1)

    with open(os.path.join(save_path, filename), 'w', newline='\n') as f:
        for line in headers:
            f.write(line)

        f.write('\n')

        diffractogram.to_csv(f, index=False, header=False, sep='\t')
            

def read_xy(data, options={}, index=0):
    
    #if 'wavelength' not in data.keys():
    # Get wavelength from scan

    if 'wavelength' in data.keys() and not type(data['wavelength']) == list:
        data['wavelength'] = [data['wavelength']]

    if not 'wavelength' in data.keys() or not data['wavelength'][index]:
        wavelength = read_metadata_from_xy(path=data['path'][index])['wavelength']
    
    else:
        wavelength = data['wavelength'][index]

    with open(data['path'][index], 'r') as f:
        position = 0

        current_line = f.readline()
        
        while current_line[0] == '#' or current_line[0] == '\'':
            position = f.tell()
            current_line = f.readline()
            
        f.seek(position)

        diffractogram = pd.read_csv(f, header=None, delim_whitespace=True)


    if diffractogram.shape[1] == 2:
        diffractogram.columns = ['2th', 'I']
    elif diffractogram.shape[1] == 3:
        diffractogram.columns = ['2th', 'I', 'sigma']


    return diffractogram, wavelength



def read_metadata_from_xy(path):

    metadata = {}
    wavelength_dict = {'Cu': 1.54059, 'Mo': 0.71073}

    with open(path, 'r') as f:
        lines = f.readlines()

        for line in lines:
            # For .xy-files output from EVA
            if 'Anode' in line:
                anode = line.split()[8].strip('"')
                metadata['wavelength'] = wavelength_dict[anode]

            
            elif 'Wavelength' in line:
                # For .xy-files output from pyFAI integration
                if line.split()[-1] == 'm':
                    metadata['wavelength'] = float(line.split()[2])*10**10

                else:
                    metadata['wavelength'] = float(line.split()[-1])


            # Get temperature - exists in .xy-files saved from HTXRD-runs in .brml-files
            if 'Temperature' in line:
                metadata['temperature'] = line.split()[-1]

            # Get timestamp - exists in .xy-files saved from .brml-files
            if 'Time' in line:
                metadata['time'] = " ".join(line.split()[2:])




    if 'wavelength' not in metadata.keys():
        metadata['wavelength'] = None
    if 'temperature' not in metadata.keys():
        metadata['temperature'] = None
    if 'time' not in metadata.keys():
        metadata['time'] = None

    return metadata


def find_wavelength_from_poni(path):

    with open(path, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if 'Wavelength' in line:
                wavelength = float(line.split()[-1])*10**10


    return wavelength




def strip_headers_from_xy(path: str, filename=None) -> None:
    ''' Strips headers from a .xy-file'''


    xy = []
    with open(path, 'r') as f:
        lines = f.readlines()

        headerlines = 0
        for line in lines:
            if line[0] == '#':
                headerlines += 1
            elif line[0] == "\'":
                headerlines += 1

            else:
                xy.append(line)

    
    if not filename:
        ext = path.split('.')[-1]
        filename = path.split(f'.{ext}')[0] + f'_noheaders.{ext}'

    
    with open(filename, 'w') as f:
        for line in xy:
            f.write(line)






def read_data(data, options={}, index=0):

    beamline_extensions = ['mar3450', 'edf', 'cbf']
    file_extension = data['path'][index].split('.')[-1]

    if file_extension in beamline_extensions:
        diffractogram, wavelength = integrate_1d(data=data, options=options, index=index)
        
    elif file_extension == 'brml':
        diffractogram, wavelength = read_brml(data=data, options=options, index=index)

    elif file_extension in['xy', 'xye']:
        diffractogram, wavelength = read_xy(data=data, options=options, index=index)

    
    if options['exclude']:

        if not isinstance(options['exclude'], list):
            options['exclude'] = [options['exclude']]

        for excl in options['exclude']:
            diffractogram['I'].loc[(diffractogram['2th'] > excl[0]) & (diffractogram['2th'] < excl[1])] = 0

   
    if options['offset'] or options['normalise']:
        # Make copy of the original intensities before any changes are made through normalisation or offset, to easily revert back if need to update.
        diffractogram['I_org'] = diffractogram['I']
        diffractogram['2th_org'] = diffractogram['2th']
        
        diffractogram = adjust_intensities(diffractogram, wavelength, index, options)


    
    diffractogram = translate_wavelengths(data=diffractogram, wavelength=wavelength)

    return diffractogram, wavelength


def adjust_intensities(diffractogram, wavelength, index, options):

    if 'current_offset_y' not in options.keys():
        options['current_offset_y'] = options['offset_y']
    else:
        if options['current_offset_y'] != options['offset_y']:
            options['offset_change'] = True
        
        options['current_offset_y'] = options['offset_y']

    options['current_offset_x'] = options['offset_x']

    

    #Apply offset along y-axis
    diffractogram['I'] = diffractogram['I_org'] # Reset intensities

    if options['normalise']:
        diffractogram['I'] = diffractogram['I'] / diffractogram['I'].max()

        if not isinstance(options['multiply'], list):
            options['multiply'] = [options['multiply']]

        diffractogram['I'] = diffractogram['I'] * options['multiply'][index]

    if options['drawdown']:
        diffractogram['I'] = diffractogram['I'] - diffractogram['I'].mean()

    diffractogram['I'] = diffractogram['I'] + index*options['offset_y']

    # Apply offset along x-axis
    relative_shift = (wavelength / 1.54059)*options['offset_x'] # Adjusts the offset-factor to account for wavelength, so that offset_x given is given in 2th_cuka-units
    diffractogram['2th'] = diffractogram['2th_org']
    diffractogram['2th'] = diffractogram['2th'] + index*relative_shift


    return diffractogram

def revert_offset(diffractogram,which=None):

    if which == 'both':
        diffractogram['2th'] = diffractogram['2th_org']
        diffractogram['I'] = diffractogram['I_org']
        
    if which == 'y':
        diffractogram['I'] = diffractogram['I_org']
    
    if which == 'x':
        diffractogram['2th'] = diffractogram['2th_org']

    return diffractogram

def load_reflection_table(data: dict, reflections_params: dict, options={}):

    required_options = ['ref_wavelength', 'to_wavelength']

    default_options = {
        'ref_wavelength': 1.54059,
        'to_wavelength': None
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    # VESTA outputs the file with a header that has a space between the parameter and units - so there is some extra code to rectify the issue
    # that ensues from this formatting
    reflections = pd.read_csv(reflections_params['path'], delim_whitespace=True)

    # Remove the extra column that appears from the headers issue
    reflections.drop(reflections.columns[-1], axis=1, inplace=True)

    with open(reflections_params['path'], 'r') as f:
        line = f.readline()

        headers = line.split()

        # Delete the fourth element which is '(Å)'
        del headers[4]

        # Change name of column to avoid using greek letters
        headers[7] = '2th'

    # Set the new modified headers as the headers of 
    reflections.columns = headers

    reflections = translate_wavelengths(data=reflections, wavelength=options['ref_wavelength'], to_wavelength=options['to_wavelength'])

    if 'heatmap' in data.keys():
        
        start_2th, stop_2th = data['diffractogram'][0]['2th'].min(), data['diffractogram'][0]['2th'].max()
        len_2th = stop_2th - start_2th 
        #print(start_2th, stop_2th, len_2th)
        
        start_heatmap, stop_heatmap = 0, data['heatmap'].shape[1]
        len_heatmap = stop_heatmap - start_heatmap
        #print(start_heatmap, stop_heatmap, len_heatmap)

        scale = len_heatmap/len_2th

        #print(scale)
        #print(stop_2th * scale)

        reflections['heatmap'] = (reflections['2th']-start_2th) * scale

    return reflections



def translate_wavelengths(data: pd.DataFrame, wavelength: float, to_wavelength=None) -> pd.DataFrame:
    # FIXME Somewhere here there is an invalid arcsin-argument. Not sure where.

    pd.options.mode.chained_assignment = None

    # Translate to CuKalpha
    cuka = 1.54059 # Å

    if cuka > wavelength:
        max_2th_cuka = 2*np.arcsin(wavelength/cuka) * 180/np.pi
    else:
        max_2th_cuka = data['2th'].max()

    data['2th_cuka'] = np.NAN

    data['2th_cuka'].loc[data['2th'] <= max_2th_cuka] = 2*np.arcsin(cuka/wavelength * np.sin((data['2th'].loc[data['2th'] <= max_2th_cuka]/2) * np.pi/180)) * 180/np.pi

    # Translate to MoKalpha
    moka = 0.71073 # Å

    if moka > wavelength:
        max_2th_moka = 2*np.arcsin(wavelength/moka) * 180/np.pi
    else:
        max_2th_moka = data['2th'].max()

    data['2th_moka'] = np.NAN

    data['2th_moka'].loc[data['2th'] <= max_2th_moka] = 2*np.arcsin(moka/wavelength * np.sin((data['2th'].loc[data['2th'] <= max_2th_moka]/2) * np.pi/180)) * 180/np.pi
    
    
    # Convert to other parameters
    data['d'] = wavelength  / (2*np.sin((2*data['2th']*np.pi/180)/2))
    data['1/d'] = 1/data['d']
    data['q'] = np.abs((4*np.pi/wavelength)*np.sin(data['2th']/2 * np.pi/180))
    data['q2'] = data['q']**2
    data['q4'] = data['q']**4


    if to_wavelength:
        if to_wavelength >= cuka:
            max_2th = 2*np.arcsin(cuka/to_wavelength) * 180/np.pi
        else:
            max_2th = data['2th_cuka'].max()

        
        data['2th'] = np.NAN
        data['2th'].loc[data['2th_cuka'] <= max_2th] = 2*np.arcsin(to_wavelength/cuka * np.sin((data['2th_cuka'].loc[data['2th_cuka'] <= max_2th]/2) * np.pi/180)) * 180/np.pi 




    return data




def trim_xy_region(path, region):

    df = pd.read_csv(path, header=None, delim_whitespace=True)
    df.columns = ['2th', 'I']

    df = df.loc[(df['2th'] > region[0]) & (df['2th'] < region[1])]
    
    folder = os.path.dirname(path)
    save_folder = os.path.join(folder, 'trimmed')

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    df.to_csv(os.path.join(save_folder, os.path.basename(path)), sep='\t', header=None, index=None)


def raise_intensities_xy(path, region=None):

    df = pd.read_csv(path, header=None, delim_whitespace=True)
    df.columns = ['2th', 'I']

    if region:
        df = df.loc[(df['2th'] > region[0]) & (df['2th'] < region[1])]

    df['I'] = df['I'] - df['I'].min()


    folder = os.path.dirname(path)
    save_folder = os.path.join(folder, 'raised')

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    df.to_csv(os.path.join(save_folder, os.path.basename(path)), sep='\t', header=None, index=None)