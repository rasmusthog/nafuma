import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def read_batsmall(path):
	''' Reads BATSMALL-data into a DataFrame.

	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	df = pd.read_csv(path, skiprows=2, sep='\t')
	df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

	return df




def read_neware(path, summary=False):
	''' Reads electrochemistry data, currently only from the Neware battery cycler. Will convert to .csv if the filetype is .xlsx,
	which is the file format the Neware provides for the backup data. In this case it matters if summary is False or not. If file
	type is .csv, it will just open the datafile and it does not matter if summary is False or not.'''
	from xlsx2csv import Xlsx2csv

	# Convert from .xlsx to .csv to make readtime faster
	if path.split('.')[-1] == 'xlsx':
		csv_details = ''.join(path.split('.')[:-1]) + '_details.csv'
		csv_summary = ''.join(path.split('.')[:-1]) + '_summary.csv'

		if not os.path.isfile(csv_summary):
			Xlsx2csv(path, outputencoding="utf-8").convert(csv_summary, sheetid=3)

		if not os.path.isfile(csv_details):
			Xlsx2csv(path, outputencoding="utf-8").convert(csv_details, sheetid=4)

		if summary:
			df = pd.read_csv(csv_summary)
		else:
			df = pd.read_csv(csv_details)

	elif path.split('.')[-1] == 'csv':
		df = pd.read_csv(path)


	return df



def read_biologic(path):
	''' Reads Bio-Logic-data into a DataFrame.
	
	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	with open(path, 'r') as f:
		lines = f.readlines()

	header_lines = int(lines[1].split()[-1]) - 1
	

	df = pd.read_csv(path, sep='\t', skiprows=header_lines)
	df.dropna(inplace=True, axis=1)

	return df







def process_batsmall_data(df, units=None, splice_cycles=None, molecular_weight=None):
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


	# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
	new_units = set_units(units=units)
	old_units = get_old_units(df, kind='batsmall')
	df = unit_conversion(df=df, new_units=new_units, old_units=old_units, kind='batsmall')

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

		sub_df = df.loc[df['count'] == i].copy()

		sub_df.loc[dchg_mask, 'current']  *= -1
		sub_df.loc[dchg_mask, 'specific_capacity'] *= -1

		chg_df = sub_df.loc[chg_mask]
		dchg_df = sub_df.loc[dchg_mask]

		# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
		if chg_df.empty and dchg_df.empty:
			continue

		cycles.append((chg_df, dchg_df))




	return cycles


def process_neware_data(df, units=None, splice_cycles=None, active_material_weight=None, molecular_weight=None, reverse_discharge=False):

	""" Takes data from NEWARE in a DataFrame as read by read_neware() and converts units, adds columns and splits into cycles.
	
	Input:
	df: pandas DataFrame containing NEWARE data as read by read_neware()
	units: dictionary containing the desired units. keywords: capacity, current, voltage, mass, energy, time
	splice_cycles: tuple containing index of cycles that should be spliced. Specifically designed to add two charge steps during the formation cycle with two different max voltages
	active_materiale_weight: weight of the active material (in mg) used in the cell. 
	molecular_weight: the molar mass (in g mol^-1) of the active material, to calculate the number of ions extracted. Assumes one electron per Li+/Na+-ion """
	
	# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
	new_units = set_units(units=units)
	old_units = get_old_units(df=df, kind='neware')
	
	df = add_columns(df=df, active_material_weight=active_material_weight, molecular_weight=molecular_weight, old_units=old_units, kind='neware')

	df = unit_conversion(df=df, new_units=new_units, old_units=old_units, kind='neware')


	# Creates masks for charge and discharge curves
	chg_mask = df['status'] == 'CC Chg'
	dchg_mask = df['status'] == 'CC DChg'

	# Initiate cycles list
	cycles = []

	# Loop through all the cycling steps, change the current and capacities in the 
	for i in range(df["cycle"].max()):

		sub_df = df.loc[df['cycle'] == i].copy()

		#sub_df.loc[dchg_mask, 'current']  *= -1
		#sub_df.loc[dchg_mask, 'capacity'] *= -1

		chg_df = sub_df.loc[chg_mask]
		dchg_df = sub_df.loc[dchg_mask]

		# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
		if chg_df.empty and dchg_df.empty:
			continue

		if reverse_discharge:
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


def process_biologic_data(df, units=None, splice_cycles=None, active_material_weight=None, molecular_weight=None, reverse_discharge=False):

	# Pick out necessary columns
	df = df[['Ns changes', 'Ns', 'time/s', 'Ewe/V', 'Energy charge/W.h', 'Energy discharge/W.h', '<I>/mA',  'Capacity/mA.h', 'cycle number']].copy()

	# Complete set of new units and get the units used in the dataset, and convert values in the DataFrame from old to new.
	new_units = set_units(units=units)
	old_units = get_old_units(df=df, kind='biologic')
	
	df = add_columns(df=df, active_material_weight=active_material_weight, molecular_weight=molecular_weight, old_units=old_units, kind='biologic')

	df = unit_conversion(df=df, new_units=new_units, old_units=old_units, kind='biologic')


	# Creates masks for charge and discharge curves
	chg_mask = (df['status'] == 1) & (df['status_change'] != 1)
	dchg_mask = (df['status'] == 2) & (df['status_change'] != 1)



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

		if reverse_discharge:
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


def add_columns(df, active_material_weight, molecular_weight, old_units, kind):

	if kind == 'neware':
		if active_material_weight:
			df["SpecificCapacity({}/mg)".format(old_units["capacity"])] = df["Capacity({})".format(old_units['capacity'])] / (active_material_weight)

			if molecular_weight:
				faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
				seconds_per_hour = 3600 # s h^-1
				f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

				df["IonsExtracted"] = (df["SpecificCapacity({}/mg)".format(old_units['capacity'])]*molecular_weight)*1000/f


	if kind == 'biologic':
		if active_material_weight:

			capacity = old_units['capacity'].split('h')[0] + '.h'


			df["SpecificCapacity({}/mg)".format(old_units["capacity"])] = df["Capacity/{}".format(capacity)] / (active_material_weight)

			if molecular_weight:
				faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
				seconds_per_hour = 3600 # s h^-1
				f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

				df["IonsExtracted"] = (df["SpecificCapacity({}/mg)".format(old_units['capacity'])]*molecular_weight)*1000/f

	return df


def unit_conversion(df, new_units, old_units, kind):
	from . import unit_tables

	if kind == 'batsmall':

		df["TT [{}]".format(old_units["time"])] = df["TT [{}]".format(old_units["time"])] * unit_tables.time()[old_units["time"]].loc[new_units['time']]
		df["U [{}]".format(old_units["voltage"])] = df["U [{}]".format(old_units["voltage"])] * unit_tables.voltage()[old_units["voltage"]].loc[new_units['voltage']]
		df["I [{}]".format(old_units["current"])] = df["I [{}]".format(old_units["current"])] * unit_tables.current()[old_units["current"]].loc[new_units['current']]
		df["C [{}/{}]".format(old_units["capacity"], old_units["mass"])] = df["C [{}/{}]".format(old_units["capacity"], old_units["mass"])] * (unit_tables.capacity()[old_units["capacity"]].loc[new_units["capacity"]] / unit_tables.mass()[old_units["mass"]].loc[new_units["mass"]])

		df.columns = ['time', 'voltage', 'current', 'count', 'specific_capacity', 'comment']


	if kind == 'neware':
		df['Current({})'.format(old_units['current'])] = df['Current({})'.format(old_units['current'])] * unit_tables.current()[old_units['current']].loc[new_units['current']]
		df['Voltage({})'.format(old_units['voltage'])] = df['Voltage({})'.format(old_units['voltage'])] * unit_tables.voltage()[old_units['voltage']].loc[new_units['voltage']]
		df['Capacity({})'.format(old_units['capacity'])] = df['Capacity({})'.format(old_units['capacity'])] * unit_tables.capacity()[old_units['capacity']].loc[new_units['capacity']]
		df['Energy({})'.format(old_units['energy'])] = df['Energy({})'.format(old_units['energy'])] * unit_tables.energy()[old_units['energy']].loc[new_units['energy']]
		df['CycleTime({})'.format(new_units['time'])] = df.apply(lambda row : convert_time_string(row['Relative Time(h:min:s.ms)'], unit=new_units['time']), axis=1)
		df['RunTime({})'.format(new_units['time'])] = df.apply(lambda row : convert_datetime_string(row['Real Time(h:min:s.ms)'], reference=df['Real Time(h:min:s.ms)'].iloc[0], unit=new_units['time']), axis=1)
		columns = ['status', 'jump', 'cycle', 'steps', 'current', 'voltage', 'capacity', 'energy']

		if 'SpecificCapacity({}/mg)'.format(old_units['capacity']) in df.columns:
			df['SpecificCapacity({}/mg)'.format(old_units['capacity'])] = df['SpecificCapacity({}/mg)'.format(old_units['capacity'])] * unit_tables.capacity()[old_units['capacity']].loc[new_units['capacity']] / unit_tables.mass()['mg'].loc[new_units["mass"]]
			columns.append('specific_capacity')

			if 'IonsExtracted' in df.columns:
				columns.append('ions')



		columns.append('cycle_time')
		columns.append('run_time')


		df.drop(['Record number', 'Relative Time(h:min:s.ms)', 'Real Time(h:min:s.ms)'], axis=1, inplace=True)

		df.columns = columns

	if kind == 'biologic':
		df['time/{}'.format(old_units['time'])] = df["time/{}".format(old_units["time"])] * unit_tables.time()[old_units["time"]].loc[new_units['time']]
		df["Ewe/{}".format(old_units["voltage"])] = df["Ewe/{}".format(old_units["voltage"])] * unit_tables.voltage()[old_units["voltage"]].loc[new_units['voltage']]
		df["<I>/{}".format(old_units["current"])] = df["<I>/{}".format(old_units["current"])] * unit_tables.current()[old_units["current"]].loc[new_units['current']]
		
		capacity = old_units['capacity'].split('h')[0] + '.h'
		df["Capacity/{}".format(capacity)] = df["Capacity/{}".format(capacity)] * (unit_tables.capacity()[old_units["capacity"]].loc[new_units["capacity"]])

		columns = ['status_change', 'status', 'time', 'voltage', 'energy_charge', 'energy_discharge', 'current', 'capacity', 'cycle']

		if 'SpecificCapacity({}/mg)'.format(old_units['capacity']) in df.columns:
			df['SpecificCapacity({}/mg)'.format(old_units['capacity'])] = df['SpecificCapacity({}/mg)'.format(old_units['capacity'])] * unit_tables.capacity()[old_units['capacity']].loc[new_units['capacity']] / unit_tables.mass()['mg'].loc[new_units["mass"]]
			columns.append('specific_capacity')

			if 'IonsExtracted' in df.columns:
				columns.append('ions')

		df.columns = columns

	return df


def set_units(units=None):
	
	# Complete the list of units - if not all are passed, then default value will be used
	required_units = ['time', 'current', 'voltage', 'capacity', 'mass', 'energy']
	default_units = {'time': 'h', 'current': 'mA', 'voltage': 'V', 'capacity': 'mAh', 'mass': 'g', 'energy': 'mWh'}

	if not units:
		units = default_units

	if units:
		for unit in required_units:
			if unit not in units.keys():
				units[unit] = default_units[unit]


	return units



def get_old_units(df, kind):
	
	if kind=='batsmall':
		time = df.columns[0].split()[-1].strip('[]')
		voltage = df.columns[1].split()[-1].strip('[]')
		current = df.columns[2].split()[-1].strip('[]')
		capacity, mass = df.columns[4].split()[-1].strip('[]').split('/')
		old_units = {'time': time, 'current': current, 'voltage': voltage, 'capacity': capacity, 'mass': mass}

	if kind=='neware':
		
		for column in df.columns:
			if 'Voltage' in column:
				voltage = column.split('(')[-1].strip(')')
			elif 'Current' in column:
				current = column.split('(')[-1].strip(')')
			elif 'Capacity' in column:
				capacity = column.split('(')[-1].strip(')')
			elif 'Energy' in column:
				energy = column.split('(')[-1].strip(')')

		old_units = {'voltage': voltage, 'current': current, 'capacity': capacity, 'energy': energy}


	if kind=='biologic':

		for column in df.columns:
			if 'time' in column:
				time = column.split('/')[-1]
			elif 'Ewe' in column:
				voltage = column.split('/')[-1]
			elif 'Capacity' in column:
				capacity = column.split('/')[-1].replace('.', '')
			elif 'Energy' in column:
				energy = column.split('/')[-1].replace('.', '')
			elif '<I>' in column:
				current = column.split('/')[-1]

		old_units = {'voltage': voltage, 'current': current, 'capacity': capacity, 'energy': energy, 'time': time}
	


	return old_units

def convert_time_string(time_string, unit='ms'):
	''' Convert time string from Neware-data with the format hh:mm:ss.xx to any given unit'''

	h, m, s = time_string.split(':')  
	ms = float(s)*1000 + int(m)*1000*60 + int(h)*1000*60*60

	factors = {'ms': 1, 's': 1/1000, 'min': 1/(1000*60), 'h': 1/(1000*60*60)}

	t = ms*factors[unit]

	return t



def convert_datetime_string(datetime_string, reference, unit='s'):
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

	time = s * factors[unit]

	return time

def splice_cycles(first, second):

	first_chg = first[0]
	first_dchg = first[1]
	first

	second_chg = second[0]
	second_dchg = second[1]

	chg_df = first[0].append(second[0])

	return True








