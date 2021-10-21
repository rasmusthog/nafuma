import fabio, pyFAI
import numpy as np
import os


def get_image_array(path):

    image = fabio.open(path)
    image_array = image.data

    return image_array


def integrate_1d(path, calibrant, bins, options):

    required_options = ['unit', 'extension']

    default_options = {'unit': '2th_deg', 'extension': '_integrated.dat'}

    if not options:
        options = default_options
    
    else:
        for option in required_options:
            if option not in options.keys():
                options[option] = default_options[option]
    
    
    image = get_image_array(path)
    ai = pyFAI.load(calibrant)

    filename = os.path.split(path)[-1].split('.')[0] + options['extension']

    res = ai.integrate1d(image, bins, unit=options['unit'], filename=filename)

    



def view_integrator(calibrant):
    ''' Prints out information about the azimuthal integrator
    
    Input:
    calibrant: Path to the azimuthal integrator file (.PONI)
    
    Output:
    None'''

    ai = pyFAI.load(calibrant)

    print("pyFAI version:", pyFAI.version)
    print("\nIntegrator: \n", ai)


