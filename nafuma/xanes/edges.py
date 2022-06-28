import pandas as pd
import numpy as np
from scipy.constants import c, h

# From 2019 redefinition of SI base units: https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units
keV_per_J = (1 / 1.602176634e-19) / 1000

# kXu values taken from International Tables for Crystallography Volume , Kulwer Academic Publishers - Dordrect / Boston / London (1992)
k_edge = {   'Z':    [ 1, 2,
                        3, 4, 5, 6, 7, 8, 9, 10,
                        11, 12, 13, 14, 15, 16, 17, 18,
                        19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
                        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48], 
            'Atom': [   'H', 'He',
                        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'], 
            'kXu': [    np.nan, np.nan,
                        226.5, np.nan, np.nan, 43.68, 30.99, 23.32, np.nan, np.nan, 
                        np.nan, 9.5117, 7.9511, 6.7446, 5.7866, 5.0182, 4.3969, 3.8707, 
                        3.43645, 3.07016, 2.7573, 2.49730, 2.26902, 2.07012, 1.89636, 1.74334, 1.60811, 1.48802, 1.38043, 1.2833, 1.19567, 1.11652, 1.04497, 0.97978, 0.91995, 0.86547,
                        0.81549, 0.76969, 0.72762, 0.68877, 0.65291, 0.61977, 0.5891, 0.56047, 0.53378, 0.50915, 0.48582, 0.46409]}


k_edge = pd.DataFrame(k_edge)
k_edge['keV'] = np.round(h*c/(k_edge['kXu']*10**-10) * keV_per_J, 3)