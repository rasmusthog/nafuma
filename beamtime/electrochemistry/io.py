import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_battsmall(path):
	''' Reads BATTSMALL-data into a DataFrame.

	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	df = pd.read_csv(path, skiprows=2, sep='\t')
	df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

	return df




def read_neware(path, summary=False, active_material_weight=None, molecular_weight=None):
    ''' Reads electrochemistry data, currently only from the Neware battery cycler. Will convert to .csv if the filetype is .xlsx,
    which is the file format the Neware provides for the backup data. In this case it matters if summary is False or not. If file
    type is .csv, it will just open the datafile and it does not matter if summary is False or not.'''
    

    # Convert from .xlsx to .csv to make readtime faster
    if filename.split('.')[-1] == 'xlsx':
        csv_details = ''.join(filename.split('.')[:-1]) + '_details.csv'
        csv_summary = ''.join(filename.split('.')[:-1]) + '_summary.csv'

        Xlsx2csv(filename, outputencoding="utf-8").convert(csv_summary, sheetid=3)
        Xlsx2csv(filename, outputencoding="utf-8").convert(csv_details, sheetid=4)
    
        if summary:
            df = pd.read_csv(csv_summary)
        else:
            df = pd.read_csv(csv_details)

    elif filename.split('.')[-1] == 'csv':

        df = pd.read_csv(filename)
    
    
    return df







#def process_battsmall_data(df, t='ms', C='mAh/g', I='mA', U='V'):

def process_battsmall_data(df, units=None, splice_cycles=None, molecular_weight=None):
	''' Takes BATTSMALL-data in the form of a DataFrame and cleans the data up and converts units into desired units.
	Splits up into individual charge and discharge DataFrames per cycle, and outputs a list where each element is a tuple with the Chg and DChg-data. E.g. cycles[10][0] gives the charge data for the 11th cycle.

	For this to work, the cycling program must be set to use the counter.

	Input:
	df (required): A pandas DataFrame containing BATTSMALL-data, as obtained from read_battsmall().
	t (optional): Unit for time data. Defaults to ms.
	C (optional): Unit for specific capacity. Defaults to mAh/g.
	I (optional): Unit for current. Defaults mA.
	U (optional): Unit for voltage. Defaults to V.

	Output:
	cycles: A list with 
	'''

	#########################
	#### UNIT CONVERSION ####
	#########################

	# Complete the list of units - if not all are passed, then default value will be used
	required_units = ['t', 'I', 'U', 'C']
	default_units = {'t': 'h', 'I': 'mA', 'U': 'V', 'C': 'mAh/g'}

	if not units:
		units = default_units

	if units:
		for unit in required_units:
			if unit not in units.values():
				units[unit] = default_units[unit]
		
	
	# Get the units used in the data set
	t_prev = df.columns[0].split()[-1].strip('[]')
	U_prev = df.columns[1].split()[-1].strip('[]')
	I_prev = df.columns[2].split()[-1].strip('[]')
	C_prev, m_prev = df.columns[4].split()[-1].strip('[]').split('/')
	prev_units = {'t': t_prev, 'I': I_prev, 'U': U_prev, 'C': C_prev}

	# Convert all units to the desired units.
	df = unit_conversion(df=df, units=units)

	# Replace NaN with empty string in the Comment-column and then remove all steps where the program changes - this is due to inconsistent values for current  
	df[["Comment"]] = df[["Comment"]].fillna(value={'Comment': ''})
	df = df[df["Comment"].str.contains("program")==False]

	# Creates masks for charge and discharge curves
	chg_mask = df['I'] >= 0
	dchg_mask = df['I'] < 0

	# Initiate cycles list
	cycles = []

	# Loop through all the cycling steps, change the current and capacities in the 
	for i in range(df["Z1"].max()):

		sub_df = df.loc[df['Z1'] == i].copy()

		sub_df.loc[dchg_mask, 'I']  *= -1
		sub_df.loc[dchg_mask, 'C'] *= -1

		chg_df = sub_df.loc[chg_mask]
		dchg_df = sub_df.loc[dchg_mask]

		# Continue to next iteration if the charge and discharge DataFrames are empty (i.e. no current)
		if chg_df.empty and dchg_df.empty:
			continue

		cycles.append((chg_df, dchg_df))




	return cycles


def process_neware_data(df, units=None, splice_cycles=None, active_material_weight=None, molecular_weight=None):

	#########################
	#### UNIT CONVERSION ####
	#########################
	
	# Complete the list of units - if not all are passed, then default value will be used
	required_units = ['t', 'I', 'U', 'C']
	default_units = {'t': 'h', 'I': 'mA', 'U': 'V', 'C': 'mAh/g'}

	if not units:
		units = default_units

	if units:
		for unit in required_units:
			if unit not in units.values():
				units[unit] = default_units[unit]

	

	# Get the units used in the data set
	t_prev = 's' # default in 
	U_prev = df.columns[1].split()[-1].strip('[]')
	I_prev = df.columns[2].split()[-1].strip('[]')
	C_prev, m_prev = df.columns[4].split()[-1].strip('[]').split('/')
	prev_units = {'t': t_prev, 'I': I_prev, 'U': U_prev, 'C': C_prev}

	# Convert all units to the desired units.
	df = unit_conversion(df=df, units=units)



	if active_material_weight:
		df["SpecificCapacity(mAh/g)"] = df["Capacity(mAh)"] / (active_material_weight / 1000)

	if molecular_weight:
		faradays_constant = 96485.3365 # [F] = C mol^-1 = As mol^-1
		seconds_per_hour = 3600 # s h^-1
		f = faradays_constant / seconds_per_hour * 1000.0 # [f] = mAh mol^-1

		df["IonsExtracted"] = (df["SpecificCapacity(mAh/g)"]*molecular_weight)/f


def unit_conversion(df, units, prev_units, kind):

	C, m = units['C'].split('/')
	C_prev, m_prev = prev_units['C'].split('/')


	# Define matrix for unit conversion for time
	t_units_df = {'h': [1, 60, 3600, 3600000], 'min': [1/60, 1, 60, 60000], 's': [1/3600, 1/60, 1, 1000], 'ms': [1/3600000, 1/60000, 1/1000, 1]}
	t_units_df = pd.DataFrame(t_units_df)
	t_units_df.index = ['h', 'min', 's', 'ms']

	# Define matrix for unit conversion for current
	I_units_df = {'A': [1, 1000, 1000000], 'mA': [1/1000, 1, 1000], 'uA': [1/1000000, 1/1000, 1]}
	I_units_df = pd.DataFrame(I_units_df)
	I_units_df.index = ['A', 'mA', 'uA']

	# Define matrix for unit conversion for voltage
	U_units_df = {'V': [1, 1000, 1000000], 'mV': [1/1000, 1, 1000], 'uV': [1/1000000, 1/1000, 1]}
	U_units_df = pd.DataFrame(U_units_df)
	U_units_df.index = ['V', 'mV', 'uV']

	# Define matrix for unit conversion for capacity
	C_units_df = {'Ah': [1, 1000, 1000000], 'mAh': [1/1000, 1, 1000], 'uAh': [1/1000000, 1/1000, 1]}
	C_units_df = pd.DataFrame(C_units_df)
	C_units_df.index = ['Ah', 'mAh', 'uAh']

	# Define matrix for unit conversion for capacity
	m_units_df = {'kg': [1, 1000, 1000000, 1000000000], 'g': [1/1000, 1, 1000, 1000000], 'mg': [1/1000000, 1/1000, 1, 1000], 'ug': [1/1000000000, 1/1000000, 1/1000, 1]}
	m_units_df = pd.DataFrame(m_units_df)
	m_units_df.index = ['kg', 'g', 'mg', 'ug']

	#print(df["TT [{}]".format(t_prev)])
	df["TT [{}]".format(t_prev)] = df["TT [{}]".format(t_prev)] * t_units_df[t_prev].loc[units['t']]
	df["U [{}]".format(U_prev)] = df["U [{}]".format(U_prev)] * U_units_df[U_prev].loc[units['U']]
	df["I [{}]".format(I_prev)] = df["I [{}]".format(I_prev)] * I_units_df[I_prev].loc[units['I']]
	df["C [{}/{}]".format(C_prev, m_prev)] = df["C [{}/{}]".format(C_prev, m_prev)] * (C_units_df[C_prev].loc[C] / m_units_df[m_prev].loc[m])

	df.columns = ['TT', 'U', 'I', 'Z1', 'C', 'Comment']



	return df



def convert_time_string(time_string, unit='ms'):
	''' Convert time string from Neware-data with the format hh:mm:ss.xx to any given unit'''

	h, m, s = time_string.split(':')  
	ms = int(s)*1000 + int(m)*1000*60 + int(h)*1000*60*60

	factors = {'ms': 1, 's': 1/1000, 'min': 1/(1000*60), 'h': 1/(1000*60*60)}

	t = ms*factors[unit]

	return t



def convert_datetime_string(datetime_string, reference, unit='s'):
	''' Convert time string from Neware-data with the format yyy-mm-dd hh:mm:ss to any given unit'''

	from datetime import datetime

	# Parse the 
	cur_date, cur_time = datetime_string.split()
	cur_y, cur_mo, cur_d = cur_date.split('-')
	cur_h, cur_m, cur_s = cur_time.split(':')
	cur_date = datetime(int(cur_y), int(cur_mo), int(cur_d), int(cur_h), int(cur_m), int(cur_s))

	ref_date, ref_time = reference.split()
	ref_y, ref_mo, ref_d = ref_date.split('-')
	ref_h, ref_m, ref_s = ref_time.split(':')
	ref_date = datetime(int(ref_y), int(ref_mo), int(ref_d), int(ref_h), int(ref_m), int(ref_s))

	days = cur_date - ref_date

	s = days.seconds

	factors = {'ms': 1000, 's': 1, 'min': 1/(60), 'h': 1/(60*60)}

	t = s * factors[unit]

	return t

def splice_cycles(first, second):

	first_chg = first[0]
	first_dchg = first[1]
	first

	second_chg = second[0]
	second_dchg = second[1]

	chg_df = first[0].append(second[0])

	return True








