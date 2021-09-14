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



def unit_conversion(df, units):

	C, m = units['C'].split('/')

	# Get the units used in the data set
	t_prev = df.columns[0].split()[-1].strip('[]')
	U_prev = df.columns[1].split()[-1].strip('[]')
	I_prev = df.columns[2].split()[-1].strip('[]')
	C_prev, m_prev = df.columns[4].split()[-1].strip('[]').split('/')


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
	df["C [{}/{}]".format(C_prev, m_prev)] = df["C [{}/{}]".format(C_prev, m_prev)] * (C_units_df[C_prev].loc[units['C']] / m_units_df[m_prev].loc[units['m']])

	df.columns = ['TT', 'U', 'I', 'Z1', 'C', 'Comment']



	return df

#def process_battsmall_data(df, t='ms', C='mAh/g', I='mA', U='V'):

def process_battsmall_data(df, units=None):
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

	required_units = ['t', 'I', 'U', 'C']
	default_units = {'t': 'h', 'I': 'mA', 'U': 'V', 'C': 'mAh/g'}

	if not units:
		units = default_units

	if units:
		for unit in required_units:
			if unit not in units.values():
				units[unit] = default_units[unit]
			

	# Convert all units to the desired units.
	df = unit_conversion(df=df, units=units)

	# Replace NaN with empty string in the Comment-column and then remove all steps where the program changes - this is due to inconsistent values for current and 
	df[["Comment"]] = df[["Comment"]].fillna(value={'Comment': ''})
	df = df[df["Comment"].str.contains("program")==False]

	# Creates masks for 
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


		cycles.append((chg_df, dchg_df))


	return cycles



def plot_gc(cycles, which_cycles='all', chg=True, dchg=True, colours=None, x='C', y='U'):

	fig, ax = prepare_gc_plot()


	if which_cycles == 'all':
		which_cycles = [i for i, c in enumerate(cycles)]

	if not colours:
		chg_colour = (40/255, 70/255, 75/255) # Dark Slate Gray #28464B
		dchg_colour = (239/255, 160/255, 11/255) # Marigold #EFA00B
	
	

	for i, cycle in cycles:
		if i in which_cycles:
			if chg:
				cycle[0].plot(ax=ax)




	






def prepare_gc_plot(figsize=(14,7), dpi=None):

	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
	

	return fig, ax






