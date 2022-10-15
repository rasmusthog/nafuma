import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import nafuma.auxillary as aux
from sympy import re


# FIXME This is not good practice, but a temporary fix as I don't have time to understand what causes the SettingWithCopyWarning.
# Read this: https://www.dataquest.io/blog/settingwithcopywarning/
pd.set_option('mode.chained_assignment', None)


def read_data(data, options={}):

	if data['kind'] == 'neware':
		df = read_neware(data['path'], options=options)
		cycles = process_neware_data(df=df, options=options)

	elif data['kind'] == 'batsmall':
		df = read_batsmall(data['path'])
		cycles = process_batsmall_data(df=df, options=options)

	elif data['kind'] == 'biologic':
		df = read_biologic(data['path'])
		cycles = process_biologic_data(df=df, options=options)

	return cycles



def read_neware(path, options={}):
	''' Reads electrochemistry data, currently only from the Neware battery cycler. Will convert to .csv if the filetype is .xlsx,
	which is the file format the Neware provides for the backup data. In this case it matters if summary is False or not. If file
	type is .csv, it will just open the datafile and it does not matter if summary is False or not.'''
	from xlsx2csv import Xlsx2csv

	# FIXME Do a check if a .csv-file already exists even if the .xlsx is passed

	# Convert from .xlsx to .csv to make readtime faster
	if path.endswith('xlsx'):
		csv_details = ''.join(path[:-5]) + '_details.csv'
		csv_summary = ''.join(path[:-5]) + '_summary.csv'

	
		if not os.path.isfile(csv_summary):
			Xlsx2csv(path, outputencoding="utf-8").convert(csv_summary, sheetid=3)

		if not os.path.isfile(csv_details):
			Xlsx2csv(path, outputencoding="utf-8").convert(csv_details, sheetid=4)

		if options['summary']:
			df = pd.read_csv(csv_summary)
		else:
			df = pd.read_csv(csv_details)

	elif path.endswith('csv'):
		df = pd.read_csv(path)

	return df


def read_batsmall(path):
	''' Reads BATSMALL-data into a DataFrame.

	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''


	# Determine if decimal point is . or ,
	with open(path, 'r') as f:
		for i, line in enumerate(f):
			if i == 10:
				values = line.split()
				if len(values[1].split('.')) == 2:
					decimal_point = '.'
				elif len(values[1].split(',')) == 2:
					decimal_point = ','
				


	df = pd.read_csv(path, skiprows=2, sep='\t', decimal=decimal_point)
	df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

	return df


def read_biologic(path):
	''' Reads Bio-Logic-data into a DataFrame.
	
	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	with open(path, 'rb') as f:
		lines = f.readlines()

	header_lines = int(lines[1].split()[-1]) - 1
	

	df = pd.read_csv(path, sep='\t', skiprows=header_lines, encoding='cp1252')
	df.dropna(inplace=True, axis=1)

	return df



def process_batsmall_data(df, options=None):
	''' Takes BATSMALL-data in the form of a DataFrame and cleans the data up and converts units into desired units.
	Splits up into individual charge and discharge DataFrames per cycle, and outputs a list where each element is a tuple with the Chg and DChg-data. E.g. cycles[10][0] gives the charge data for the 11th cycle.

	For this to work, the cycling program must be set to use the counter.

	Input:
	df (required): A pandas DataFrame containing BATSMALL-data, as obtained from read_batsmall().
	t (optional): Unit for time data. Defaults to ms.
	C (optional): Unit for specific capacity. Defaults to mAh/g.
	I (optional): Unit for current. Defaults mA.
	U (optional): Unit for voltage. Defaults to V.

	Output:
	cycles: A list with 
	'''

	required_options = ['splice_cycles', 'molecular_weight', 'reverse_discharge', 'units']
	
	default_options = {
		'splice_cycles': False, 
		'molecular_weight': None, 
		'reverse_discharge': False, 
		'units': None}


	aux.update_options(options=options, required_options=required_options, default_options=default_options)
	options['kind'] = 'batsmall'

	# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
	set_units(options)
	options['old_units'] = get_old_units(df, options)


	df = add_columns(df=df, options=options) # adds columns to the DataFrame if active material weight and/or molecular weight has been passed in options

	df = unit_conversion(df=df, options=options)


	if options['splice_cycles']:
		df = splice_cycles(df=df, options=options)


	# Replace NaN with empty string in the Comment-column and then remove all steps where the program changes - this is due to inconsistent values for current  
	df[["comment"]] = df[["comment"]].fillna(value={'comment': ''})
	df = df[df["comment"].str.contains("program")==False]

	# Creates masks for charge and discharge curves
	chg_mask = df['current'] >= 0
	dchg_mask = df['current'] < 0

	# Initiate cycles list
	cycles = []

	# Loop through all the cycling steps, change the current and capacities in the 
	for i in range(df["count"].max()):

		sub_df = df.loc[df['count'] == i+1].copy()

		sub_df.loc[dchg_mask, 'current']  *= -1
		sub_df.loc[dchg_mask, 'specific_capacity'] *= -1

		chg_df = sub_df.loc[chg_mask]
		dchg_df = sub_df.loc[dchg_mask]

		# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
		if chg_df.empty and dchg_df.empty:
			continue

		chg_df['reaction_coordinate'] = chg_df['time'] * np.abs(chg_df['current'].mean())
		dchg_df['reaction_coordinate'] = dchg_df['time'] * np.abs(dchg_df['current'].mean())

		if options['reverse_discharge']:
			max_capacity = dchg_df['capacity'].max() 
			dchg_df['capacity'] = np.abs(dchg_df['capacity'] - max_capacity)

			if 'specific_capacity' in df.columns:
				max_capacity = dchg_df['specific_capacity'].max() 
				dchg_df['specific_capacity'] = np.abs(dchg_df['specific_capacity'] - max_capacity)

				if 'ions' in df.columns:
					max_capacity = dchg_df['ions'].max() 
					dchg_df['ions'] = np.abs(dchg_df['ions'] - max_capacity)

		cycles.append((chg_df, dchg_df))


	return cycles


def splice_cycles(df, options: dict) -> pd.DataFrame:
	''' Splices two cycles together - if e.g. one charge cycle are split into several cycles due to change in parameters. 
	
	Incomplete, only accomodates BatSmall so far, and only for charge.'''

	if options['kind'] == 'batsmall':

		# Creates masks for charge and discharge curves
		chg_mask = df['current'] >= 0

		# Loop through all the cycling steps, change the current and capacities in the 
		for i in range(df["count"].max()):
			sub_df = df.loc[df['count'] == i+1]
			sub_df_chg = sub_df.loc[chg_mask]

			# get indices where the program changed
			chg_indices = sub_df_chg[sub_df_chg["comment"].str.contains("program")==True].index.to_list()
			
			# Delete first item if first cycle after rest (this will just be the start of the cycling)
			if i+1 == 1:
				del chg_indices[0]
			
	
			if chg_indices:
				last_chg =  chg_indices.pop()
					
			if chg_indices:
				for i in chg_indices:
					add = df['specific_capacity'].iloc[i-1]
					df['specific_capacity'].iloc[i:last_chg] = df['specific_capacity'].iloc[i:last_chg] + add


	if options['kind'] == 'neware':
		
		if options['summary']:
			for i in range(df['cycle'].max()):
				sub_df = df.loc[df['cycle'] == i+1].copy()

				if sub_df['status'].loc[sub_df['status'] == 'CC Chg'].count() > 1:
					indices = sub_df.index[sub_df['status'] == 'CC Chg']
					
					add_columns = ['capacity', 'specific_capacity', 'ions', 'energy', 'cycle_time']
					
					for column in add_columns:
						if column in df.columns:
							df[column].iloc[indices[-1]] = df[column].iloc[indices[-1]] + df[column].iloc[indices[0]]

					df.drop(index=indices[0], inplace=True)
					df.reset_index(inplace=True, drop=True)

		else:
			for i in range(df['cycle'].max()):
				sub_df = df.loc[df['cycle'] == i+1].copy()
				sub_chg_df = sub_df.loc[sub_df['status'] == 'CC Chg'].copy()

				steps_indices = sub_chg_df['steps'].unique()

				if len(steps_indices) > 1:

					add_columns = ['capacity', 'specific_capacity', 'ions', 'energy', 'cycle_time']

					for column in add_columns:
						if column in df.columns:
							# Extract the maximum value from the first of the two cycles by accessing the column value of the highest index of the first cycle
							add = df[column].iloc[df.loc[df['steps'] == steps_indices[0]].index.max()]
							
							df[column].loc[df['steps'] == steps_indices[1]] += add

	return df




def process_neware_data(df, options={}):

	""" Takes data from NEWARE in a DataFrame as read by read_neware() and converts units, adds columns and splits into cycles.
	
	Input:
	df: pandas DataFrame containing NEWARE data as read by read_neware()
	units: dictionary containing the desired units. keywords: capacity, current, voltage, mass, energy, time
	splice_cycles: tuple containing index of cycles that should be spliced. Specifically designed to add two charge steps during the formation cycle with two different max voltages
	active_materiale_weight: weight of the active material (in mg) used in the cell. 
	molecular_weight: the molar mass (in g mol^-1) of the active material, to calculate the number of ions extracted. Assumes one electron per Li+/Na+-ion """
	
	required_options = ['units', 'active_material_weight', 'molecular_weight', 'reverse_discharge', 'splice_cycles', 'increment_cycles_from', 'delete_datapoints']
	
	default_options = {
		'units': None, 
		'active_material_weight': None, 
		'molecular_weight': None, 
		'reverse_discharge': False, 
		'splice_cycles': None,
		'increment_cycles_from': None,# index
		'delete_datapoints': None, # list of indices
		}  


	aux.update_options(options=options, required_options=required_options, default_options=default_options)
	options['kind'] = 'neware'


	if not options['summary']:
		# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
		set_units(options=options) # sets options['units']
		options['old_units'] = get_old_units(df=df, options=options)
		
		df = add_columns(df=df, options=options) # adds columns to the DataFrame if active material weight and/or molecular weight has been passed in options

		df = unit_conversion(df=df, options=options) # converts all units from the old units to the desired units

		if options['increment_cycles_from']:
			df['cycle'].iloc[options['increment_cycles_from']:] += 1

		if options['delete_datapoints']:
			for datapoint in options['delete_datapoints']:
				df.drop(index=datapoint, inplace=True)

		if options['splice_cycles']:
			df = splice_cycles(df=df, options=options)


		# Creates masks for charge and discharge curves
		chg_mask = df['status'] == 'CC Chg'
		dchg_mask = df['status'] == 'CC DChg'

		# Initiate cycles list
		cycles = []



		# Loop through all the cycling steps, change the current and capacities in the 
		for i in range(df["cycle"].max()):

			sub_df = df.loc[df['cycle'] == i+1].copy()

			chg_df = sub_df.loc[chg_mask]
			dchg_df = sub_df.loc[dchg_mask]

			# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
			if chg_df.empty and dchg_df.empty:
				continue


			# Reverses the discharge curve if specified
			if options['reverse_discharge']:
				max_capacity = dchg_df['capacity'].max() 
				dchg_df['capacity'] = np.abs(dchg_df['capacity'] - max_capacity)

				if 'specific_capacity' in df.columns:
					max_capacity = dchg_df['specific_capacity'].max() 
					dchg_df['specific_capacity'] = np.abs(dchg_df['specific_capacity'] - max_capacity)

					if 'ions' in df.columns:
						max_capacity = dchg_df['ions'].max() 
						dchg_df['ions'] = np.abs(dchg_df['ions'] - max_capacity)


			if not chg_df.empty:
				chg_df.reset_index(inplace=True)
			if not dchg_df.empty:
				dchg_df.reset_index(inplace=True)
			
			cycles.append((chg_df, dchg_df))


		return cycles


	elif options['summary']:
		set_units(options=options)
		options['old_units'] = get_old_units(df=df, options=options)

		df = add_columns(df=df, options=options)
		df = unit_conversion(df=df, options=options)
		

		if options['splice_cycles']:
			df = splice_cycles(df=df, options=options)


		chg_df = df.loc[df['status'] == 'CC Chg']
		chg_df.reset_index(inplace=True)
		dchg_df = df.loc[df['status'] == 'CC DChg']
		dchg_df.reset_index(inplace=True)

		# Construct new DataFrame
		new_df = pd.DataFrame(chg_df["cycle"])
		new_df.insert(1,'charge_capacity',chg_df['capacity']) 
		new_df.insert(1,'charge_specific_capacity',chg_df['specific_capacity'])
		new_df.insert(1,'discharge_capacity',dchg_df['capacity'])
		new_df.insert(1,'discharge_specific_capacity',dchg_df['specific_capacity'])
		new_df.insert(1,'charge_energy',chg_df['energy'])
		new_df.insert(1,'charge_specific_energy',chg_df['specific_energy'])
		new_df.insert(1,'discharge_energy',dchg_df['energy'])
		new_df.insert(1,'discharge_specific_energy',dchg_df['specific_energy'])

		new_df = calculate_efficiency(df=new_df, options=options)

		return new_df



def process_biologic_data(df, options=None):

	required_options = ['units', 'active_material_weight', 'molecular_weight', 'reverse_discharge', 'splice_cycles']
	
	default_options = {
		'units': None, 
		'active_material_weight': None, 
		'molecular_weight': None, 
		'reverse_discharge': False, 
		'splice_cycles': None}


	# Check if the DataFrame contains GC or CV data.
	# FIXME This might not be a very rigorous method of checking. E.g. Rest has mode == 3, so if loading a short GC with many Rest-datapoints, the mean will be 2 and it will be treated as CV. For now manual override is sufficient
	if not 'mode' in options.keys():
		options['mode'] = 'GC' if int(df['mode'].mean()) == 1 else 'CV'
		
	aux.update_options(options=options, required_options=required_options, default_options=default_options)
	options['kind'] = 'biologic'

	# Pick out necessary columns
	headers = [
		'Ns changes', 'Ns', 'time/s', 'Ewe/V', 'Energy charge/W.h', 'Energy discharge/W.h', '<I>/mA',  'Capacity/mA.h', 'cycle number' ] if options['mode'] == 'GC' else [
		'ox/red', 'time/s', 'control/V', 'Ewe/V', '<I>/mA', 'cycle number', '(Q-Qo)/C', 'P/W'
		]


	df = df[headers].copy()

	# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
	set_units(options)
	options['old_units'] = get_old_units(df=df, options=options)
	
	df = add_columns(df=df, options=options)

	df = unit_conversion(df=df, options=options)


	# Creates masks for charge and discharge curves
	if options['mode'] == 'GC':
		chg_mask = (df['status'] == 1) & (df['status_change'] != 1)
		dchg_mask = (df['status'] == 2) & (df['status_change'] != 1)

	elif options['mode'] == 'CV':
		chg_mask = (df['status'] == 1) # oxidation
		dchg_mask = (df['status'] == 0) # reduction

	# Initiate cycles list
	cycles = []

	# Loop through all the cycling steps, change the current and capacities in the 
	for i in range(int(df["cycle"].max())):

		sub_df = df.loc[df['cycle'] == i].copy()

		#sub_df.loc[dchg_mask, 'current']  *= -1
		#sub_df.loc[dchg_mask, 'capacity'] *= -1

		chg_df = sub_df.loc[chg_mask]
		dchg_df = sub_df.loc[dchg_mask]

		# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
		if chg_df.empty and dchg_df.empty:
			continue

		if options['mode'] == 'GC' and options['reverse_discharge']:
			max_capacity = dchg_df['capacity'].max() 
			dchg_df['capacity'] = np.abs(dchg_df['capacity'] - max_capacity)

			if 'specific_capacity' in df.columns:
				max_capacity = dchg_df['specific_capacity'].max() 
				dchg_df['specific_capacity'] = np.abs(dchg_df['specific_capacity'] - max_capacity)

				if 'ions' in df.columns:
					max_capacity = dchg_df['ions'].max() 
					dchg_df['ions'] = np.abs(dchg_df['ions'] - max_capacity)

		
		if options['mode'] == 'CV':
			chg_df = chg_df.sort_values(by='voltage').reset_index(drop=True)
			dchg_df = dchg_df.sort_values(by='voltage', ascending=False).reset_index(drop=True)
		
		cycles.append((chg_df, dchg_df))


	return cycles


def add_columns(df, options):
	from . import unit_tables

	if options['kind'] == 'neware':

		if options['summary']:
			df[f'Energy({options["old_units"]["energy"]})'] = np.abs(df[f'Net discharge energy({options["old_units"]["energy"]})'])

		if options['active_material_weight']:
			df[f"SpecificCapacity({options['old_units']['capacity']}/mg)"] = df["Capacity({})".format(options['old_units']['capacity'])] / (options['active_material_weight'])
			df[f"SpecificEnergy({options['old_units']['energy']}/mg)"] = df["Energy({})".format(options['old_units']['energy'])] / (options['active_material_weight'])

			if options['molecular_weight']:
				faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
				seconds_per_hour = 3600 # s h^-1
				f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

				df["IonsExtracted"] = (df["SpecificCapacity({}/mg)".format(options['old_units']['capacity'])]*options['molecular_weight'])*1000/f


	if options['kind'] == 'biologic':
		if options['active_material_weight']:

			capacity = options['old_units']['capacity'].split('h')[0] + '.h'


			df["SpecificCapacity({}/mg)".format(options['old_units']["capacity"])] = df["Capacity/{}".format(capacity)] / (options['active_material_weight'])

			if options['molecular_weight']:
				faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
				seconds_per_hour = 3600 # s h^-1
				f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

				df["IonsExtracted"] = (df["SpecificCapacity({}/mg)".format(options['old_units']['capacity'])]*options['molecular_weight'])*1000/f


	if options['kind'] == 'batsmall':
		if options['active_material_weight']:


			active_material_weight = options['active_material_weight'] * unit_tables.mass()['mg'].loc[options['units']['mass']]
			capacity = options['old_units']['capacity']

			df[f'Capacity [{options["old_units"]["capacity"]}]'] = df[f'C [{options["old_units"]["capacity"]}/{options["old_units"]["mass"]}]'] * active_material_weight

			if options['molecular_weight']:
				faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
				seconds_per_hour = 3600 # s h^-1
				f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

				molecular_weight = options['molecular_weight'] * unit_tables.mass()['g'].loc[options['units']['mass']]
				df["IonsExtracted"] = (df[f'C [{options["old_units"]["capacity"]}/{options["old_units"]["mass"]}]'] * molecular_weight)/f


		#df['reaction_coordinate'] = (df[f'TT [{options["old_units"]["time"]}]'] * unit_tables.time()[options['old_units']["time"]].loc["h"]) / np.abs(df[f'I [{options["old_units"]["current"]}]'] * unit_tables.current()[options['old_units']["current"]].loc['A'])


	return df


def calculate_efficiency(df: pd.DataFrame, options: dict) -> pd.DataFrame:

	
	default_options = {
		'reference_index': 0
	}

	options = aux.update_options(options=options, required_options=default_options.keys(), default_options=default_options)

	df['charge_capacity_fade'] 		= (df['charge_capacity'] / df['charge_capacity'].iloc[options['reference_index']])*100
	df['discharge_capacity_fade'] 	= (df['discharge_capacity'] / df['discharge_capacity'].iloc[options['reference_index']])*100

	df['coulombic_efficiency'] 		= (df['discharge_capacity'] / df['charge_capacity'])*100
	df['energy_efficiency'] 		= (df['discharge_energy'] / df['charge_energy'])*100


	return df


def unit_conversion(df, options):
	from . import unit_tables

	if options['kind'] == 'batsmall':

		df["TT [{}]".format(options['old_units']["time"])] = df["TT [{}]".format(options['old_units']["time"])] * unit_tables.time()[options['old_units']["time"]].loc[options['units']['time']]
		df["U [{}]".format(options['old_units']["voltage"])] = df["U [{}]".format(options['old_units']["voltage"])] * unit_tables.voltage()[options['old_units']["voltage"]].loc[options['units']['voltage']]
		df["I [{}]".format(options['old_units']["current"])] = df["I [{}]".format(options['old_units']["current"])] * unit_tables.current()[options['old_units']["current"]].loc[options['units']['current']]
		df["C [{}/{}]".format(options['old_units']["capacity"], options['old_units']["mass"])] = df["C [{}/{}]".format(options['old_units']["capacity"], options['old_units']["mass"])] * (unit_tables.capacity()[options['old_units']["capacity"]].loc[options['units']["capacity"]] / unit_tables.mass()[options['old_units']["mass"]].loc[options['units']["mass"]])

		columns = ['time', 'voltage', 'current', 'count', 'specific_capacity', 'comment']

		# Add column labels for specific capacity and ions if they exist
		if f'Capacity [{options["old_units"]["capacity"]}]' in df.columns:
			df[f'Capacity [{options["old_units"]["capacity"]}]'] = df[f'Capacity [{options["old_units"]["capacity"]}]'] * unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']]
			columns.append('capacity')

		if 'IonsExtracted' in df.columns:
			columns.append('ions')

		#columns.append('reaction_coordinate')
		
		df.columns = columns





	if options['kind'] == 'neware':


		record_number = 'Data serial number' if 'Data serial number' in df.columns else 'Record number'
		relative_time = 'Relative Time(h:min:s.ms)' if 'Relative Time(h:min:s.ms)' in df.columns else 'Relative Time'
		continuous_time = 'Continuous Time(h:min:s.ms)' if 'Continuous Time(h:min:s.ms)' in df.columns else 'Continuous Time'
		real_time = 'Real Time(h:min:s.ms)' if 'Real Time(h:min:s.ms)' in df.columns else 'Real Time'
	
		
		if options['summary']:
			df[f'Energy({options["old_units"]["energy"]})'] 			= df[f'Energy({options["old_units"]["energy"]})']				* unit_tables.energy()[options['old_units']['energy']].loc[options['units']['energy']]
			df[f'Starting current({options["old_units"]["current"]})'] 	= df[f'Starting current({options["old_units"]["current"]})'] 	* unit_tables.current()[options['old_units']['current']].loc[options['units']['current']]
			df[f'Start Volt({options["old_units"]["voltage"]})'] 		= df[f'Start Volt({options["old_units"]["voltage"]})'] 			* unit_tables.voltage()[options['old_units']['voltage']].loc[options['units']['voltage']]
			df[f'Capacity({options["old_units"]["capacity"]})'] 		= df[f'Capacity({options["old_units"]["capacity"]})'] 			* unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']]
			df[f'Energy({options["old_units"]["energy"]})'] 			= df[f'Energy({options["old_units"]["energy"]})'] 				* unit_tables.energy()[options['old_units']['energy']].loc[options['units']['energy']]
			df[f'CycleTime({options["units"]["time"]})'] 				= df.apply(lambda row : convert_time_string(row[relative_time], unit=options['units']['time']), axis=1)
			df[f'RunTime({options["units"]["time"]})'] 					= df.apply(lambda row : convert_datetime_string(row[real_time], reference=df[real_time].iloc[0], ref_time=df[f'CycleTime({options["units"]["time"]})'].iloc[0],unit=options['units']['time']), axis=1)



			droplist = [
				'Chnl',
				'Original step',
				f'End Volt({options["old_units"]["voltage"]})',
				f'Termination current({options["old_units"]["current"]})',
				relative_time,
				real_time,
				continuous_time,
				f'Net discharge capacity({options["old_units"]["capacity"]})',
				f'Chg Cap({options["old_units"]["capacity"]})',
				f'DChg Cap({options["old_units"]["capacity"]})',
				f'Net discharge energy({options["old_units"]["energy"]})',
				f'Chg Eng({options["old_units"]["energy"]})', 
				f'DChg Eng({options["old_units"]["energy"]})'
			]

			# Drop all undesireable columns
			for drop in droplist:
				if drop in df.columns:
					df.drop(drop, axis=1, inplace=True)

			columns = ['cycle', 'steps', 'status', 'voltage', 'current', 'capacity', 'energy']


			# Add column labels for specific capacity and ions if they exist
			if 'SpecificCapacity({}/mg)'.format(options['old_units']['capacity']) in df.columns:
				df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] = df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] * unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']] / unit_tables.mass()['mg'].loc[options['units']["mass"]]
				columns.append('specific_capacity')

			if f'SpecificEnergy({options["old_units"]["energy"]}/mg)' in df.columns:
				df[f'SpecificEnergy({options["old_units"]["energy"]}/mg)'] = df[f'SpecificEnergy({options["old_units"]["energy"]}/mg)'] * unit_tables.energy()[options['old_units']['energy']].loc[options['units']['energy']] / unit_tables.mass()['mg'].loc[options['units']["mass"]]
				columns.append('specific_energy')

			if 'IonsExtracted' in df.columns:
				columns.append('ions')

			# Append energy column label here as it was the last column to be generated
			columns.append('cycle_time')
			columns.append('runtime')

			# Apply new column labels 
			df.columns = columns
	

		else:
			df['Current({})'.format(options['old_units']['current'])] = df['Current({})'.format(options['old_units']['current'])] * unit_tables.current()[options['old_units']['current']].loc[options['units']['current']]
			df['Voltage({})'.format(options['old_units']['voltage'])] = df['Voltage({})'.format(options['old_units']['voltage'])] * unit_tables.voltage()[options['old_units']['voltage']].loc[options['units']['voltage']]
			df['Capacity({})'.format(options['old_units']['capacity'])] = df['Capacity({})'.format(options['old_units']['capacity'])] * unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']]
			df['Energy({})'.format(options['old_units']['energy'])] = df['Energy({})'.format(options['old_units']['energy'])] * unit_tables.energy()[options['old_units']['energy']].loc[options['units']['energy']]
			df['CycleTime({})'.format(options['units']['time'])] = df.apply(lambda row : convert_time_string(row[relative_time], unit=options['units']['time']), axis=1)
			df['RunTime({})'.format(options['units']['time'])] = df.apply(lambda row : convert_datetime_string(row[real_time], reference=df[real_time].iloc[0],  ref_time=df[f'CycleTime({options["units"]["time"]})'].iloc[0], unit=options['units']['time']), axis=1)
			columns = ['status', 'jump', 'cycle', 'steps', 'current', 'voltage', 'capacity', 'energy']

			if 'SpecificCapacity({}/mg)'.format(options['old_units']['capacity']) in df.columns:
				df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] = df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] * unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']] / unit_tables.mass()['mg'].loc[options['units']["mass"]]
				columns.append('specific_capacity')

			if f'SpecificEnergy({options["old_units"]["energy"]}/mg)' in df.columns:
				df[f'SpecificEnergy({options["old_units"]["energy"]}/mg)'] = df[f'SpecificEnergy({options["old_units"]["energy"]}/mg)'] * unit_tables.energy()[options['old_units']['energy']].loc[options['units']['energy']] / unit_tables.mass()['mg'].loc[options['units']["mass"]]
				columns.append('specific_energy')


			if 'IonsExtracted' in df.columns:
				columns.append('ions')

			columns.append('cycle_time')
			columns.append('time')


			droplist = [record_number, relative_time, real_time]

			for drop in droplist:
				if drop in df.columns:
					df.drop(drop, axis=1, inplace=True)

			df.columns = columns

	if options['kind'] == 'biologic':
		for column in df.columns:
			if 'time' in column:
				df['time/{}'.format(options['old_units']['time'])] = df["time/{}".format(options['old_units']["time"])] * unit_tables.time()[options['old_units']["time"]].loc[options['units']['time']]
			
			if 'Ewe' in column:
				df["Ewe/{}".format(options['old_units']["voltage"])] = df["Ewe/{}".format(options['old_units']["voltage"])] * unit_tables.voltage()[options['old_units']["voltage"]].loc[options['units']['voltage']]

			if '<I>' in column:
				df["<I>/{}".format(options['old_units']["current"])] = df["<I>/{}".format(options['old_units']["current"])] * unit_tables.current()[options['old_units']["current"]].loc[options['units']['current']]
		
			if 'Capacity' in column:
				capacity = options['old_units']['capacity'].split('h')[0] + '.h'
				df["Capacity/{}".format(capacity)] = df["Capacity/{}".format(capacity)] * (unit_tables.capacity()[options['old_units']["capacity"]].loc[options['units']["capacity"]])

		

		columns = [
			'status_change', 'status', 'time', 'voltage', 'energy_charge', 'energy_discharge', 'current', 'capacity', 'cycle'] if options['mode'] == 'GC' else [	# GC headers
			'status', 'time', 'control_voltage', 'voltage', 'current', 'cycle', 'charge', 'power' 	# CV headers
		]

		if options['mode'] == 'GC':
			if 'SpecificCapacity({}/mg)'.format(options['old_units']['capacity']) in df.columns:
				df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] = df['SpecificCapacity({}/mg)'.format(options['old_units']['capacity'])] * unit_tables.capacity()[options['old_units']['capacity']].loc[options['units']['capacity']] / unit_tables.mass()['mg'].loc[options['units']["mass"]]
				columns.append('specific_capacity')

				if 'IonsExtracted' in df.columns:
					columns.append('ions')

		df.columns = columns

	return df


def set_units(options: dict) -> None:
	
	# Complete the list of units - if not all are passed, then default value will be used
	required_units = ['time', 'current', 'voltage', 'capacity', 'mass', 'energy', 'specific_capacity']
	
	default_units = {
		'time': 'h', 
		'current': 'mA', 
		'voltage': 'V', 
		'capacity': 'mAh', 
		'mass': 'g', 
		'energy': 'mWh', 
		'specific_capacity': None}

	if not options['units']:
		options['units'] = default_units


	aux.update_options(options=options['units'], required_options=required_units, default_options=default_units)

	options['units']['specific_capacity'] = r'{} {}'.format(options['units']['capacity'], options['units']['mass']) + '$^{-1}$'



def get_old_units(df: pd.DataFrame, options: dict) -> dict:
	''' Reads a DataFrame with cycling data and determines which units have been used and returns these in a dictionary'''
	
	if options['kind'] == 'batsmall':

		old_units = {}

		for column in df.columns:
			if 'TT [' in column:
				old_units['time'] = column.split()[-1].strip('[]')
			elif 'U [' in column:
				old_units['voltage'] = column.split()[-1].strip('[]')
			elif 'I [' in column:
				old_units['current'] = column.split()[-1].strip('[]')
			elif 'C [' in column:
				old_units['capacity'], old_units['mass'] = column.split()[-1].strip('[]').split('/')

		# time = df.columns[0].split()[-1].strip('[]')
		# voltage = df.columns[1].split()[-1].strip('[]')
		# current = df.columns[2].split()[-1].strip('[]')
		# capacity, mass = df.columns[4].split()[-1].strip('[]').split('/')
		# old_units = {'time': time, 'current': current, 'voltage': voltage, 'capacity': capacity, 'mass': mass}

	if options['kind']=='neware':

		for column in df.columns:
			if 'Voltage' in column or 'Start Volt' in column:
				voltage = column.split('(')[-1].strip(')')
			elif 'Current' in column or 'Starting current' in column:
				current = column.split('(')[-1].strip(')')
			elif 'Capacity' in column:
				capacity = column.split('(')[-1].strip(')')
			elif 'Energy' in column or 'Eng' in column:
				energy = column.split('(')[-1].strip(')')

		old_units = {'voltage': voltage, 'current': current, 'capacity': capacity, 'energy': energy}


	if options['kind'] == 'biologic':

		old_units = {}
		for column in df.columns:
			if 'time' in column:
				old_units['time'] = column.split('/')[-1]
			elif 'Ewe' in column:
				old_units['voltage'] = column.split('/')[-1]
			elif 'Capacity' in column:
				old_units['capacity'] = column.split('/')[-1].replace('.', '')
			elif 'Energy' in column:
				old_units['energy'] = column.split('/')[-1].replace('.', '')
			elif '<I>' in column:
				old_units['current'] = column.split('/')[-1]
	
	return old_units

def convert_time_string(time_string, unit='ms'):
	''' Convert time string from Neware-data with the format hh:mm:ss.xx to any given unit'''

	h, m, s = time_string.split(':')
	ms = float(s)*1000 + int(m)*1000*60 + int(h)*1000*60*60

	factors = {'ms': 1, 's': 1/1000, 'min': 1/(1000*60), 'h': 1/(1000*60*60)}

	t = ms*factors[unit]

	return t



def convert_datetime_string(datetime_string, reference, ref_time, unit='s'):
	''' Convert time string from Neware-data with the format yyy-mm-dd hh:mm:ss to any given unit'''

	from datetime import datetime

	# Parse the 
	current_date, current_time = datetime_string.split()
	current_year, current_month, current_day = current_date.split('-')
	current_hour, current_minute, current_second = current_time.split(':')
	current_date = datetime(int(current_year), int(current_month), int(current_day), int(current_hour), int(current_minute), int(current_second))

	reference_date, reference_time = reference.split()
	reference_year, reference_month, reference_day = reference_date.split('-')
	reference_hour, reference_minute, reference_second = reference_time.split(':')
	reference_date = datetime(int(reference_year), int(reference_month), int(reference_day), int(reference_hour), int(reference_minute), int(reference_second))

	days = current_date - reference_date


	s = days.days*24*60*60 + days.seconds

	factors = {'ms': 1000, 's': 1, 'min': 1/(60), 'h': 1/(60*60)}

	time = s * factors[unit] + ref_time

	return time






