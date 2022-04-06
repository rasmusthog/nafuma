import beamtime.plotting as btp
from cycler import cycler
import itertools
import numpy as np

import matplotlib.pyplot as plt


def test_generate_colours() -> None:

    assert type(btp.generate_colours('black', kind='single')) == itertools.cycle

    palettes = [('qualitative', 'Dark2_8')]
    colour_cycle = btp.generate_colours(palettes)

    assert type(colour_cycle) == itertools.cycle


    # Test that it actually loaded 8 colours when given a set of 8 colours to 

    same_colour = None
    for i in range(10):
        colour = next(colour_cycle)
        if i == 0:
            first_colour = colour
        
        if colour == first_colour:
            repeat_colour_index = i

        
    assert repeat_colour_index == 8




def test_update_rc_params() -> None:

    rc_params = {
        'lines.linewidth': 100
    }

    prev_params = plt.rcParams['lines.linewidth']
    
    # Update run commands if any is passed (will pass an empty dictionary if not passed)
    btp.update_rc_params(rc_params)

    new_params = plt.rcParams['lines.linewidth']
    
    assert new_params == 100
    assert prev_params != new_params


    # Reset run commands
    plt.rcdefaults()



def test_scale_figure() -> None:

    width, height = 1, 1

    format_params = {
        'upscaling_factor': 2,
        'compress_width': 1,
        'compress_height': 1
    }

    width1, height1 = btp.scale_figure(format_params=format_params, width=width, height=height)

    assert width1 == 2 and height1 == 2

    format_params = {
        'upscaling_factor': 1,
        'compress_width': 0.5,
        'compress_height': 1
    }

    width2, height2 = btp.scale_figure(format_params=format_params, width=width, height=height)

    assert width2 == 0.5 and height2 == 1

    format_params = {
        'upscaling_factor': 2,
        'compress_width': 0.5,
        'compress_height': 0.2
    }

    width2, height2 = btp.scale_figure(format_params=format_params, width=width, height=height)

    assert width2 == 1 and height2 == 0.4

    
def test_determine_width() -> None:

    conversion_cm_inch = 0.3937008 # cm to inch

    format_params = {
        'column_type': 'single',
        'single_column_width': 5,
        'double_column_width': 10,
        'width_ratio': '1:1'
    }

    assert np.round(btp.determine_width(format_params),6) == np.round(5*conversion_cm_inch,6)
    
    format_params['column_type'] = 'double'

    assert np.round(btp.determine_width(format_params), 6) == np.round(10*conversion_cm_inch, 6)
    

    format_params['column_type'] = 'single'
    format_params['width_ratio'] = '1:2'

    assert np.round(btp.determine_width(format_params), 6) == np.round(2.5*conversion_cm_inch, 6)

def test_determine_height() -> None:


    width = 1

    format_params = {
        'aspect_ratio': '1:1'
    }

    assert btp.determine_height(format_params=format_params, width=width) == 1

    format_params['aspect_ratio'] = '3:1'

    assert (btp.determine_height(format_params=format_params, width=width) - 0.333333333333333) < 10e-7

    assert True

