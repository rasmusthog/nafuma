from pickle import MARK
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

import pandas as pd
import numpy as np
import math
import os
import shutil
from PIL import Image

import ipywidgets as widgets
from IPython.display import display

import nafuma.electrochemistry as ec
import nafuma.plotting as btp
import nafuma.auxillary	as aux




def plot_gc(data, options=None):


	# Update options
	required_options = ['force_reload', 'x_vals', 'y_vals', 'which_cycles', 'limit', 'exclude_cycles', 'show_plot', 'summary', 'charge', 'discharge', 'colours', 'differentiate_charge_discharge', 'gradient', 'interactive', 'interactive_session_active', 'rc_params', 'format_params', 'save_gif', 'save_path', 'fps']	
	default_options = {
		'force_reload': False,
		'x_vals': 'capacity', 'y_vals': 'voltage', 
		'which_cycles': 'all',
		'limit': None, # Limit line to be drawn
		'exclude_cycles': [],
		'show_plot': True,
		'summary': False,
		'charge': True, 'discharge': True, 
		'colours': None, 
		'differentiate_charge_discharge': True, 
		'gradient': False,
		'interactive': False,
		'interactive_session_active': False, 
		'rc_params': {},
		'format_params': {},
		'save_gif': False,
		'save_path': 'animation.gif',
		'fps': 1
		}

	options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

	
	# Read data if not already loaded
	if not 'cycles' in data.keys() or options['force_reload']:
		data['cycles'] = ec.io.read_data(data=data, options=options)

	
	# Update list of cycles to correct indices
	update_cycles_list(data=data, options=options)

	if options['interactive']:
		options['interactive'], options['interactive_session_active'] = False, True
		plot_gc_interactive(data=data, options=options)
		return



	colours = generate_colours(cycles=data['cycles'], options=options)
	if not options['summary']:
		
		if options['show_plot']:
			# Prepare plot
			
			if not options['fig'] and not options['ax']:
				fig, ax = btp.prepare_plot(options=options)
			else:
				fig, ax = options['fig'], options['ax']
			
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

		
		
		if options['save_gif'] and not options['interactive_session_active']:
			if not os.path.isdir('tmp'):
				os.makedirs('tmp')

			# Scale image to make GIF smaller
			options['format_params']['width'] = 7.5
			options['format_params']['height'] = 3

			options['format_params']['dpi'] = 200

			for i, cycle in enumerate(data['cycles']):
				if i in options['which_cycles']:
					
					giffig, gifax = btp.prepare_plot(options=options)

					if options['charge']:
						cycle[0].plot(x=options['x_vals'], y=options['y_vals'], ax=gifax, c=colours[i][0])
					if options['discharge']:
						cycle[1].plot(x=options['x_vals'], y=options['y_vals'], ax=gifax, c=colours[i][1])

					gifax.text(x=gifax.get_xlim()[1]*0.8, y=3, s=f'{i+1}')
					update_labels(options)

					giffig, gifax = btp.adjust_plot(fig=giffig, ax=gifax, options=options)

					plt.savefig(os.path.join('tmp', str(i+1).zfill(4)+'.png'))
					plt.close()

				
			img_paths = [os.path.join('tmp', path) for path in os.listdir('tmp') if path.endswith('png')]
			frames = []
			for path in img_paths:
				frame = Image.open(path)
				frames.append(frame)

			frames[0].save(options['save_path'], format='GIF', append_images=frames[1:], save_all=True, duration=(1/options['fps'])*1000, loop=0)

			shutil.rmtree('tmp')




	elif options['summary'] and options['show_plot']:
		# Prepare plot
		fig, ax = btp.prepare_plot(options=options)
		
		mask = []
		for i in range(data['cycles'].shape[0]):
			if i+1 in options['which_cycles']:
				mask.append(True)
			else:
				mask.append(False)


		# Drop the last row if it is midway through a charge in order to avoid mismatch of length of mask and dataset.
		if len(mask) > data['cycles'].shape[0]:
			del mask[-1]
			data['cycles'].drop(data['cycles'].tail(1).index, inplace=True)


		# FIXME To begin, the default is that y-values correspond to x-values. This should probably be implemented in more logical and consistent manner in the future.
		if options['x_vals'] in ['coulombic_efficiency', 'energy_efficiency']:
			data['cycles'].loc[mask].plot(x='cycle', y=options['x_vals'], ax=ax, color=colours[0][0], kind='scatter', marker="$\u25EF$", s=plt.rcParams['lines.markersize'])
			if options['limit']:
				ax.axhline(y=options['limit'], ls='--', c='black')

		else:
			if options['charge']:
				yval = 'charge_' + options['x_vals']
				data['cycles'].loc[mask].plot(x='cycle', y=yval, ax=ax, color=colours[0][0], kind='scatter', marker="$\u25EF$", s=plt.rcParams['lines.markersize'])
			
			if options['discharge']:
				yval = 'discharge_' + options['x_vals']
				data['cycles'].loc[mask].plot(x='cycle', y=yval, ax=ax, color=colours[0][1], kind='scatter', marker="$\u25EF$", s=plt.rcParams['lines.markersize'])


			if options['limit']:
				ax.axhline(y=options['limit'], ls='--', c='black')


		if options['interactive_session_active']:
			update_labels(options, force=True)
		else:
			update_labels(options)

	

	if options['show_plot']:
		fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)
		return data['cycles'], fig, ax
	else:
		return data['cycles'], None, None


def plot_gc_interactive(data, options):

	w = widgets.interactive(btp.ipywidgets_update, func=widgets.fixed(plot_gc), data=widgets.fixed(data), options=widgets.fixed(options),
	charge=widgets.ToggleButton(value=True),
	discharge=widgets.ToggleButton(value=True),
	x_vals=widgets.Dropdown(options=['specific_capacity', 'capacity', 'ions', 'voltage', 'time', 'energy'], value='specific_capacity', description='X-values')
	)

	options['widget'] = w

	display(w)




def plot_cv(data, options):

	# Update options
	required_options = ['force_reload', 'x_vals', 'y_vals', 'which_cycles', 'limit', 'exclude_cycles', 'show_plot', 'charge', 'discharge', 'colours', 'differentiate_charge_discharge', 'gradient', 'interactive', 'interactive_session_active', 'rc_params', 'format_params', 'save_gif', 'save_path', 'fps']	
	default_options = {
		'force_reload': False,
		'x_vals': 'voltage', 'y_vals': 'current', 
		'which_cycles': 'all',
		'limit': None, # Limit line to be drawn
		'exclude_cycles': [],
		'show_plot': True,
		'charge': True, 'discharge': True, 
		'colours': None, 
		'differentiate_charge_discharge': True, 
		'gradient': False,
		'interactive': False,
		'interactive_session_active': False, 
		'rc_params': {},
		'format_params': {},
		'save_gif': False,
		'save_path': 'animation.gif',
		'fps': 1
		}

	options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

	
	# Read data if not already loaded
	if not 'cycles' in data.keys() or options['force_reload']:
		data['cycles'] = ec.io.read_data(data=data, options=options)

	
	# Update list of cycles to correct indices
	update_cycles_list(data=data, options=options)

	colours = generate_colours(cycles=data['cycles'], options=options)

	if options['show_plot']:
		# Prepare plot
		fig, ax = btp.prepare_plot(options=options)
		for i, cycle in enumerate(data['cycles']):
			if i in options['which_cycles']:
				if options['charge']:
					cycle[0].plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=colours[i][0])

				if options['discharge']:
					cycle[1].plot(x=options['x_vals'], y=options['y_vals'], ax=ax, c=colours[i][1])

		update_labels(options)

	
	
	if options['save_gif'] and not options['interactive_session_active']:
		if not os.path.isdir('tmp'):
			os.makedirs('tmp')

		# Scale image to make GIF smaller
		options['format_params']['width'] = 7.5
		options['format_params']['height'] = 3

		options['format_params']['dpi'] = 200

		for i, cycle in enumerate(data['cycles']):
			if i in options['which_cycles']:
				
				giffig, gifax = btp.prepare_plot(options=options)

				if options['charge']:
					cycle[0].plot(x=options['x_vals'], y=options['y_vals'], ax=gifax, c=colours[i][0])
				if options['discharge']:
					cycle[1].plot(x=options['x_vals'], y=options['y_vals'], ax=gifax, c=colours[i][1])

				gifax.text(x=gifax.get_xlim()[1]*0.8, y=3, s=f'{i+1}')
				update_labels(options)

				giffig, gifax = btp.adjust_plot(fig=giffig, ax=gifax, options=options)

				plt.savefig(os.path.join('tmp', str(i+1).zfill(4)+'.png'))
				plt.close()

			
		img_paths = [os.path.join('tmp', path) for path in os.listdir('tmp') if path.endswith('png')]
		frames = []
		for path in img_paths:
			frame = Image.open(path)
			frames.append(frame)

		frames[0].save(options['save_path'], format='GIF', append_images=frames[1:], save_all=True, duration=(1/options['fps'])*1000, loop=0)

		shutil.rmtree('tmp')



	if options['show_plot']:
		fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)
		return data['cycles'], fig, ax
	else:
		return data['cycles'], None, None

def update_labels(options, force=False):

	if 'xlabel' not in options.keys() or force:
		options['xlabel'] = options['x_vals'].capitalize().replace('_', ' ')

	if 'ylabel' not in options.keys() or force:
		options['ylabel'] = options['y_vals'].capitalize().replace('_', ' ')
		

	if 'xunit' not in options.keys() or force:
		if options['x_vals'] == 'capacity':
			options['xunit'] = options['units']['capacity']
		elif options['x_vals'] == 'specific_capacity':
			options['xunit'] = f"{options['units']['capacity']} {options['units']['mass']}$^{{-1}}$"
		elif options['x_vals'] == 'time':
			options['xunit'] = options['units']['time']
		elif options['x_vals'] == 'ions':
			options['xunit'] = None
		

	if 'yunit' not in options.keys() or force:
		if options['y_vals'] == 'voltage':
			options['yunit'] = options['units']['voltage']



		
	
def update_cycles_list(data, options: dict) -> None:

	if options['which_cycles'] == 'all':
		options['which_cycles'] = [i for i in range(len(data['cycles']))]


	elif isinstance(options['which_cycles'], list):

		cycles =[]

		for cycle in options['which_cycles']:
			if isinstance(cycle, int):
				cycles.append(cycle-1)

			elif isinstance(cycle, tuple):
				interval = [i-1 for i in range(cycle[0], cycle[1]+1)]
				cycles.extend(interval)

		
		options['which_cycles'] = cycles
	
	
	# Tuple is used to define an interval - as elements tuples can't be assigned, I convert it to a list here.
	elif isinstance(options['which_cycles'], tuple):
		which_cycles = list(options['which_cycles'])

		if which_cycles[0] <= 0:
			which_cycles[0] = 1

		elif which_cycles[1] < 0:
			which_cycles[1] = len(options['which_cycles'])


		options['which_cycles'] = [i-1 for i in range(which_cycles[0], which_cycles[1]+1)]

	
	for i, cycle in enumerate(options['which_cycles']):
		if cycle in options['exclude_cycles']:
			del options['which_cycles'][i]




def prettify_gc_plot(fig, ax, options=None):


	##################################################################
	######################### UPDATE OPTIONS #########################
	##################################################################

	# Define the required options
	required_options = [
	'columns', 
	'xticks', 'yticks', 
	'show_major_ticks', 'show_minor_ticks',
	'xlim', 'ylim',
	'hide_x_axis', 'hide_y_axis', 
	'positions',
	'x_vals', 'y_vals', 
	'xlabel', 'ylabel', 
	'units', 'sizes',
	'title'
	]


	# Define the default options
	default_options = {
		'columns': 1,
		'xticks': None, 'yticks': None,
		'show_major_ticks': [True, True, True, True], 'show_minor_ticks': [True, True, True, True],
		'xlim': None,'ylim': None,
		'hide_x_axis': False, 'hide_y_axis': False,
		'positions': {'xaxis': 'bottom', 'yaxis': 'left'},
		'x_vals': 'specific_capacity', 'y_vals': 'voltage',
		'xlabel': None, 'ylabel': None,
		'units': {'capacity': 'mAh', 'specific_capacity': r'mAh g$^{-1}$', 'time': 's', 'current': 'mA', 'energy': 'mWh', 'mass': 'g', 'voltage': 'V'},
		'sizes': None,
		'title': None 	
	}

	aux.update_options(options, required_options, default_options)


	##################################################################
	########################## DEFINE SIZES ##########################
	##################################################################

	# Define the required sizes
	required_sizes = [
		'labels', 
		'legend', 
		'title', 
		'line', 'axes', 
		'tick_labels', 
		'major_ticks', 'minor_ticks']
	
	
	
	# Define default sizes
	default_sizes = {
			'labels': 30*options['columns'],
			'legend': 30*options['columns'],
			'title': 30*options['columns'],
			'line': 3*options['columns'],
			'axes': 3*options['columns'],
			'tick_labels': 30*options['columns'],
			'major_ticks': 20*options['columns'],
			'minor_ticks': 10*options['columns']	
	}
	
	# Initialise dictionary if it doesn't exist
	if not options['sizes']:
		options['sizes'] = {}


	# Update dictionary with default values where none is supplied
	for size in required_sizes:
		if size not in options['sizes']:
			options['sizes'][size] = default_sizes[size]


	##################################################################
	########################## AXIS LABELS ###########################
	##################################################################


	if not options['xlabel']:
		options['xlabel'] = prettify_labels(options['x_vals']) + ' [{}]'.format(options['units'][options['x_vals']])
	
	else:
		options['xlabel'] = options['xlabel'] + ' [{}]'.format(options['units'][options['x_vals']])


	if not options['ylabel']:
		options['ylabel'] = prettify_labels(options['y_vals']) + ' [{}]'.format(options['units'][options['y_vals']])
	
	else:
		options['ylabel'] = options['ylabel'] + ' [{}]'.format(options['units'][options['y_vals']])

	ax.set_xlabel(options['xlabel'], size=options['sizes']['labels'])
	ax.set_ylabel(options['ylabel'], size=options['sizes']['labels'])

	##################################################################
	###################### TICK MARKS & LABELS #######################
	##################################################################

	ax.tick_params(direction='in', which='major', bottom=options['show_major_ticks'][0], left=options['show_major_ticks'][1], top=options['show_major_ticks'][2], right=options['show_major_ticks'][0], length=options['sizes']['major_ticks'], width=options['sizes']['axes'])
	ax.tick_params(direction='in', which='minor', bottom=options['show_minor_ticks'][0], left=options['show_minor_ticks'][1], top=options['show_minor_ticks'][2], right=options['show_minor_ticks'][0], length=options['sizes']['minor_ticks'], width=options['sizes']['axes'])



	# DEFINE AND SET TICK DISTANCES

	from . import unit_tables

	# Define default ticks and scale to desired units
	default_ticks = {
		'specific_capacity': [100 * (unit_tables.capacity()['mAh'].loc[options['units']['capacity']] / unit_tables.mass()['g'].loc[options['units']['mass']]), 50 * (unit_tables.capacity()['mAh'].loc[options['units']['capacity']] / unit_tables.mass()['g'].loc[options['units']['mass']])], 
		'capacity': [0.1 * (unit_tables.capacity()['mAh'].loc[options['units']['capacity']]), 0.05 * (unit_tables.capacity()['mAh'].loc[options['units']['capacity']])],
		'voltage': [0.5 * (unit_tables.voltage()['V'].loc[options['units']['voltage']]), 0.25 * (unit_tables.voltage()['V'].loc[options['units']['voltage']])],
		'time': [10 * (unit_tables.time()['h'].loc[options['units']['time']]), 5 * (unit_tables.time()['h'].loc[options['units']['time']])]
		}


	if options['positions']['yaxis'] == 'right':
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()


	# Set default tick distances for x-axis if not specified
	if not options['xticks']:

		major_xtick = default_ticks[options['x_vals']][0]
		minor_xtick = default_ticks[options['x_vals']][1]

	# Otherwise apply user input
	else:
		major_xtick = options['xticks'][0]
		minor_xtick = options['xticks'][1]


	# Set default tick distances for x-axis if not specified
	if not options['yticks']:

		major_ytick = default_ticks[options['y_vals']][0]
		minor_ytick = default_ticks[options['y_vals']][1]

	# Otherwise apply user input
	else:
		major_ytick = options['yticks'][0]
		minor_ytick = options['yticks'][1]


	# Apply values
	ax.xaxis.set_major_locator(MultipleLocator(major_xtick))
	ax.xaxis.set_minor_locator(MultipleLocator(minor_xtick))

	ax.yaxis.set_major_locator(MultipleLocator(major_ytick))
	ax.yaxis.set_minor_locator(MultipleLocator(minor_ytick))




	# SET FONTSIZE OF TICK LABELS

	plt.xticks(fontsize=options['sizes']['tick_labels'])
	plt.yticks(fontsize=options['sizes']['tick_labels'])

	##################################################################
	########################## AXES LIMITS ###########################
	##################################################################

	if options['xlim']:
		plt.xlim(options['xlim'])

	if options['ylim']:
		plt.ylim(options['ylim'])

	##################################################################
	############################# TITLE ##############################
	##################################################################

	if options['title']:
		ax.set_title(options['title'], size=options['sizes']['title'])

	##################################################################
	############################# LEGEND #############################
	##################################################################

	if ax.get_legend():
		ax.get_legend().remove()

	return fig, ax

	
def prettify_labels(label):
	
	labels_dict = {
		'capacity': 'Capacity',
		'specific_capacity': 'Specific capacity',
		'voltage': 'Voltage',
		'current': 'Current',
		'energy': 'Energy',
		'time': 'Time'
	}

	return labels_dict[label]








def generate_colours(cycles, options):

	# Assign colours from the options dictionary if it is defined, otherwise use standard colours.
	if options['colours']:
		charge_colour = options['colours'][0]
		discharge_colour = options['colours'][1]
	
	else:
		charge_colour = (40/255, 70/255, 75/255) # Dark Slate Gray #28464B, coolors.co
		discharge_colour = (239/255, 160/255, 11/255) # Marigold #EFA00B, coolors.co

		if not options['differentiate_charge_discharge']:
			discharge_colour = charge_colour



	# If gradient is enabled, find start and end points for each colour
	if options['gradient']:

		add_charge = min([(1-x)*0.75 for x in charge_colour])
		add_discharge = min([(1-x)*0.75 for x in discharge_colour])

		charge_colour_start = charge_colour
		charge_colour_end = [x+add_charge for x in charge_colour]

		discharge_colour_start = discharge_colour
		discharge_colour_end = [x+add_discharge for x in discharge_colour]



	# Generate lists of colours
	colours = []
	
	for cycle_number in range(0, len(cycles)):
		if options['gradient']:
			weight_start = (len(cycles) - cycle_number)/len(cycles)
			weight_end = cycle_number/len(cycles)

			charge_colour = [weight_start*start_colour + weight_end*end_colour for start_colour, end_colour in zip(charge_colour_start, charge_colour_end)]
			discharge_colour = [weight_start*start_colour + weight_end*end_colour for start_colour, end_colour in zip(discharge_colour_start, discharge_colour_end)]

		colours.append([charge_colour, discharge_colour])

	return colours
