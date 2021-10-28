import fabio, pyFAI
import pandas as pd
import numpy as np
import os
import shutil

import zipfile
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


    required_options = ['extract_folder']
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
    



            


