import fabio, pyFAI
import pandas as pd
import numpy as np
import os
import shutil

import zipfile
import xml.etree.ElementTree as ET


import beamtime.auxillary as aux


def get_image_array(path):

    image = fabio.open(path)
    image_array = image.data

    return image_array


def get_image_headers(path):

    image = fabio.open(path)

    return image.header


def integrate_1d(data, options={}):
    ''' Integrates an image file to a 1D diffractogram. 

    Required content of data:
    calibrant (str): path to .poni-file
    nbins (int): Number of bins to divide image into
    path (str) (optional, dependent on image): path to image file - either this or image must be specified. If both is passed, image is prioritsed
    image (NumPy 2D Array) (optional, dependent on path): image array as extracted from get_image_array

    Output:
    df: DataFrame contianing 1D diffractogram if option 'return' is True
    ''' 

    required_options = ['unit', 'save', 'save_filename', 'save_extension', 'save_folder', 'overwrite']

    default_options = {
        'unit': '2th_deg', 
        'save': False,
        'save_filename': None,
        'save_extension': '_integrated.xy',
        'save_folder': '.',
        'overwrite': False}

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    

    # Get image array from filename if not passed
    if 'image' not in data.keys():
        data['image'] = get_image_array(data['path'])
    
    # Instanciate the azimuthal integrator from pyFAI from the calibrant (.poni-file)
    ai = pyFAI.load(data['calibrant'])

    # Determine filename
    filename = make_filename(data=data, options=options)
    
    # Make save_folder if this does not exist already
    if not os.path.isdir(options['save_folder']):
        os.makedirs(options['save_folder'])


    res = ai.integrate1d(data['image'], data['nbins'], unit=options['unit'], filename=filename)

    diffractogram = read_xy(filename)

    if not options['save']:
        os.remove(filename)
        shutil.rmtree('tmp')
    
    return diffractogram
    


def make_filename(data, options):

    # Define save location for integrated diffractogram data
    if not options['save']:
            options['save_folder'] = 'tmp'
            filename = os.path.join(options['save_folder'], 'tmp_diff.dat')

    elif options['save']:

        # Case 1: No filename is given. 
        if not options['save_filename']:
            # If a path is given instead of an image array, the path is taken as the trunk of the savename
            if data['path']:
                # Make filename by joining the save_folder, the filename (with extension deleted) and adding the save_extension
                filename = os.path.join(options['save_folder'], os.path.split(data['path'])[-1].split('.')[0] + options['save_extension'])
            else:
                # Make filename just "integrated.dat" in the save_folder
                filename = os.path.join(options['save_folder'], 'integrated.xy')


        else:
            filename = os.path.join(options['save_folder'], options['save_filename'])


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


def average_images(images):
    ''' Takes a list of path to image files, reads them and averages them before returning the average image'''

    image_arrays = []

    for image in images:
        image_array = xrd.io.get_image_array(image)
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




def read_brml(path, options=None):


    required_options = ['extract_folder', 'save_folder']
    default_options = {
        'extract_folder': 'temp',
        'save_folder': None
    }


    if not options:
        options = default_options

    else:
        for option in required_options:
            if option not in options.keys():
                options[option] = default_options[option]



    if not os.path.isdir(options['extract_folder']):
        os.mkdir(options['extract_folder'])


    # Extract the RawData0.xml file from the brml-file
    with zipfile.ZipFile(path, 'r') as brml:
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

        for scantype in chain.findall('ScanInformation/ScanMode'):
            if scantype.text == 'StillScan':

                if chain.get('Description') == 'Originally measured data.':
                    for data in chain.findall('Datum'):
                        data = data.text.split(',')
                        data = [float(i) for i in data]
                        twotheta, intensity = float(data[2]), float(data[3])

        
            else:
                if chain.get('Description') == 'Originally measured data.':
                    for data in chain.findall('Datum'):
                        data = data.text.split(',')
                        twotheta, intensity = float(data[2]), float(data[3])
                        
                        if twotheta > 0:
                            diffractogram.append({'2th': twotheta, 'I': intensity})

    diffractogram = pd.DataFrame(diffractogram)



    if options['save_folder']:
        if not os.path.isdir(options['save_folder']):
            os.makedirs(options['save_folder'])

        diffractogram.to_csv(options['save_folder'])



    return diffractogram
    

def read_xy(data):

    with open(data['path'], 'r') as f:
        position = 0
    
        current_line = f.readline()
        
        while current_line[0] == '#' or "\'":
            find_wavelength(line=current_line, data=data)
            position = f.tell()
            current_line = f.readline()
            
        f.seek(position)

        diffractogram = pd.read_csv(f, header=None, delim_whitespace=True)

    if diffractogram.shape[1] == 2:
        diffractogram.columns = ['2th', 'I']
    elif diffractogram.shape[1] == 3:
        diffractogram.columns = ['2th', 'I', 'sigma']


    return diffractogram


def read_data(data, options={}):

    beamline_extensions = ['mar3450', 'edf', 'cbf']
    file_extension = data['path'].split('.')[-1]

    if file_extension in beamline_extensions:
        diffractogram = integrate_1d(data=data, options=options)
        
    elif file_extension == 'brml':
        diffractogram = read_brml(path=data['path'], options=options)

    elif file_extension in['xy', 'xye']:
        diffractogram = read_xy(data['path'])

    return diffractogram
                



def load_reflection_table(path):

    # VESTA outputs the file with a header that has a space between the parameter and units - so there is some extra code to rectify the issue
    # that ensues from this formatting
    reflections = pd.read_csv(path, delim_whitespace=True)

    # Remove the extra column that appears from the headers issue
    reflections.drop(reflections.columns[-1], axis=1, inplace=True)

    with open(path, 'r') as f:
        line = f.readline()

        headers = line.split()

        # Delete the fourth element which is '(Å)'
        del headers[4]

        # Change name of column to avoid using greek letters
        headers[7] = '2th'

    # Set the new modified headers as the headers of 
    reflections.columns = headers

    return reflections



def translate_wavelengths(diffractogram, wavelength):

    # Translate to CuKalpha
    cuka = 1.54059 # Å

    if cuka > wavelength:
        max_2th_cuka = 2*np.arcsin(wavelength/cuka) * 180/np.pi
    else:
        max_2th_cuka = diffractogram['2th'].max()

    diffractogram['2th_cuka'] = np.NAN

    diffractogram['2th_cuka'].loc[diffractogram['2th'] <= max_2th_cuka] = 2*np.arcsin(cuka/wavelength * np.sin((diffractogram['2th']/2) * np.pi/180)) * 180/np.pi

    # Translate to MoKalpha
    moka = 0.71073 # Å

    if moka > wavelength:
        max_2th_moka = 2*np.arcsin(wavelength/moka) * 180/np.pi
    else:
        max_2th_moka = diffractogram['2th'].max()

    diffractogram['2th_moka'] = np.NAN

    diffractogram['2th_moka'].loc[diffractogram['2th'] <= max_2th_moka] = 2*np.arcsin(moka/wavelength * np.sin((diffractogram['2th']/2) * np.pi/180)) * 180/np.pi
    
    
    # Convert to other parameters
    diffractogram['d'] = wavelength  / (2*np.sin(2*diffractogram['2th']) * 180/np.pi)
    diffractogram['1/d'] = 1/diffractogram['d']
    diffractogram['q'] = np.abs((4*np.pi/wavelength)*np.sin(diffractogram['2th']/2 * np.pi/180))
    diffractogram['q2'] = diffractogram['q']**2
    diffractogram['q4'] = diffractogram['q']**4



def find_wavelength(line, data):

    # Find from EVA-exports
    wavelength_dict = {'Cu': 1.54059, 'Mo': 0.71073}

    if 'Anode' in line:
        