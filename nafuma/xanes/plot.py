import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math
import datetime

#import ipywidgets as widgets
#from IPython.display import display

import nafuma.xanes as xas
import nafuma.plotting as btp
import nafuma.auxillary	as aux


def plot_xanes(data, options={}):


	# Update options
	default_options = {
		'which_scans': 'all',  # Use real numbers, not indices - update_scans_list() will adjust.
		'highlight': [],
        'xlabel': 'Energy', 'ylabel': 'Intensity',
        'xunit': 'keV', 'yunit': 'arb. u.',
        'exclude_scans': [],
		'colours': None, 
		'gradient': False,
		'rc_params': {},
		'format_params': {}}

	options = aux.update_options(options=options, default_options=default_options)

	
	if not 'xanes_data' in data.keys():
		data['xanes_data'] = xas.io.load_data(data=data, options=options)

	# Update list of cycles to correct indices
	update_scans_list(data=data, options=options)

	colours = generate_colours(scans=options['which_scans'], options=options)

	# Prepare plot, and read and process data
	
	fig, ax = btp.prepare_plot(options=options)


	# Add counter to pick out correct colour
	counter = 0
	for i, path in enumerate(data['path']):
		if i in options['which_scans']:
			lw = plt.rcParams['lines.linewidth']*5 if i in options['highlight'] else plt.rcParams['lines.linewidth']

			data['xanes_data'].plot(x='ZapEnergy', y=path, ax=ax, c=colours[counter], lw=lw)
			counter += 1


	fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)

	#if options['interactive_session_active']:
	

	return fig, ax 


def pick_out_scans(metadata: dict, timestamp: list):

	# If either start or end are None, set to way back when or way into the future
	split_operator=timestamp[0][-3] #Adding this to enable reading of both "." and ":" as operators to split hour:minute:second

	if not timestamp[0]:
		timestamp[0] = datetime.datetime.strptime('1970 01 01 00:00:00', '%Y %m %d %H:%M:%S')
	else:
		if split_operator == ".":
			timestamp[0] = datetime.datetime.strptime(timestamp[0], "%d.%b %y %H.%M.%S")
		if split_operator == ":":
			timestamp[0] = datetime.datetime.strptime(timestamp[0], "%d.%b %y %H:%M:%S")
	if not timestamp[1]:
		timestamp[1] = datetime.datetime.strptime('3000 01 01 00:00:00', '%Y %m %d %H:%M:%S')
	else:
		if split_operator == ".":
			timestamp[1] = datetime.datetime.strptime(timestamp[1], "%d.%b %y %H.%M.%S")
		if split_operator == ":":
			timestamp[1] = datetime.datetime.strptime(timestamp[1], "%d.%b %y %H:%M:%S")

	scans = []
	for i, time in enumerate(metadata['time']):
		if time >= timestamp[0] and time <= timestamp[1]:
			scans.append(i)


	return scans

		




def update_scans_list(data, options: dict) -> None:

	if options['which_scans'] == 'all':
		options['which_scans'] = [i for i in range(len(data['path']))]


	elif isinstance(options['which_scans'], list):

		scans =[]

		for scan in options['which_scans']:
			if isinstance(scan, int):
				scans.append(scan-1)

			elif isinstance(scan, tuple):
				interval = [i-1 for i in range(scan[0], scan[1]+1)]
				scans.extend(interval)

		
		options['which_scans'] = scans
	
	
	# Tuple is used to define an interval - as elements tuples can't be assigned, I convert it to a list here.
	elif isinstance(options['which_scans'], tuple):
		which_scans = list(options['which_scans'])

		if which_scans[0] <= 0:
			which_scans[0] = 1

		elif which_scans[1] < 0:
			which_scans[1] = len(options['which_scans'])


		options['which_scans'] = [i-1 for i in range(which_scans[0], which_scans[1]+1)]

	
	for i, scan in enumerate(options['which_scans']):
		if scan in options['exclude_scans']:
			del options['which_scans'][i]




def generate_colours(scans, options):
    # FIXME Make this a generalised function and use this instead of this and in the electrochemsitry submodule

	# Assign colours from the options dictionary if it is defined, otherwise use standard colours.
	if options['colours']:
		colour = options['colours']
	
	else:
		#colour =  (214/255, 143/255, 214/255) # Plum Web (#D68FD6), coolors.co
		colour = (90/255, 42/255, 39/255) # Caput Mortuum(#5A2A27), coolors.co

	# If gradient is enabled, find start and end points for each colour
	if options['gradient']:

		if isinstance(colour, list) and len(colour) == 2:
			options['number_of_colours'] = len(scans)
			colours = btp.mix_colours(colour1=colour[0], colour2=colour[1], options=options)


		else:
			add = min([(1-x)*0.75 for x in colour])

			colour_start = colour
			colour_end = [x+add for x in colour]
			


	# Generate lists of colours
	if not isinstance(colour, list):
		colours = []
		for scan_number in range(0, len(scans)):

			if options['gradient']:
				weight_start = (len(scans) - scan_number)/len(scans)
				weight_end = scan_number/len(scans)

				colour = [weight_start*start_colour + weight_end*end_colour for start_colour, end_colour in zip(colour_start, colour_end)]
			
			colours.append(colour)
			
	return colours