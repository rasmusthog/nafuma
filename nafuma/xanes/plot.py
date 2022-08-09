import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math

#import ipywidgets as widgets
#from IPython.display import display

import nafuma.xanes as xas
import nafuma.plotting as btp
import nafuma.auxillary	as aux


def plot_xanes(data, options=None):


	# Update options
	required_options = ['which_scans', 'colours', 'gradient', 'rc_params', 'format_params']	
	default_options = {
		'which_scans': 'all', 
		'colours': None, 
		'gradient': False,
		'rc_params': {},
		'format_params': {}}

	options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

	
	if not 'xanes_data' in data.keys():
		data['xanes_data'] = xas.io.load_data(data=data, options=options)

	# Update list of cycles to correct indices
	update_scans_list(scans=data['path'], options=options)

	colours = generate_colours(cycles=data['cycles'], options=options)

	if options['interactive']:
		options['interactive'], options['interactive_session_active'] = False, True
		plot_gc_interactive(data=data, options=options)
		return


	# Prepare plot, and read and process data
	
	fig, ax = btp.prepare_plot(options=options)

	for i, cycle in enumerate(data['cycles']):
		if i in options['which_cycles']:
			if options['charge']:
				cycle[0].plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=colours[i][0])

			if options['discharge']:
				cycle[1].plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=colours[i][1])


	if options['interactive_session_active']:
		update_labels(options, force=True)
	else:
		update_labels(options)

	fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)

	#if options['interactive_session_active']:
	

	return data['cycles'], fig, ax 





def update_scans_list(scans, options: dict) -> None:

	if options['which_scans'] == 'all':
		options['which_scans'] = [i for i in range(len(scans))]


	elif type(options['which_scans']) == list:
		options['which_scans'] = [i-1 for i in options['which_scans']]
	
	
	# Tuple is used to define an interval - as elements tuples can't be assigned, I convert it to a list here.
	elif type(options['which_cycles']) == tuple:
		which_cycles = list(options['which_cycles'])

		if which_cycles[0] <= 0:
			which_cycles[0] = 1

		elif which_cycles[1] < 0:
			which_cycles[1] = len(cycles)


		options['which_cycles'] = [i-1 for i in range(which_cycles[0], which_cycles[1]+1)]