import pandas as pd


def read_battsmall(path):
	''' Reads BATTSMALL-data into a DataFrame.

	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	df = pd.read_csv(data1, skiprows=2, sep='\t')
	df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

	return df


def clean_battsmall_data(df, t='ms', C='mAh/g', I='mA', U='V'):
	''' Takes BATTSMALL-data in the form of a DataFrame and cleans the data up and converts units into desired units.
	Also adds a column indicating whether or not it is charge or discharge.

	Input:
	df (required): A pandas DataFrame containing BATTSMALL-data, as obtained from read_battsmall().
	t (optional): Unit for time data. Defaults to ms.
	C (optional): Unit for specific capacity. Defaults to mAh/g.
	I (optional): Unit for current. Defaults mA.
	U (optional): Unit for voltage. Defaults to V.

	Output:
	df: A cleaned up DataFrame with desired units and additional column about charge/discharge.
	'''








def unit_conversion(df, t, C, I, U):


	# Get the units used in the data set
	t_prev = df.columns[0].split()[-1].strip('[]')
	U_prev = df.columns[1].split()[-1].strip('[]')
	I_prev = df.columns[2].split()[-1].strip('[]')
	C_prev_num, C_prev_den = df.columns[4].split()[-1].strip('[]').split('/')

	# Define matrix for unit conversion for time
	t_units_df = {'h': [1, 60, 3600, 3600000], 'min': [1/60, 1, 60, 60000], 's': [1/3600, 1/60, 1, 1000], 'ms': [1/3600000, 1/60000, 1/1000, 1]}
	t_units_df = pd.DataFrame(t_units_df)
	t_units_df.index = ['h', 'min', 's', 'ms']

	# Define matrix for unit conversion for current
	I_units_df = {'A': [1, 60, 3600, 3600000], 'mA': [1/60, 1, 60, 60000], 'uA': [1/3600, 1/60, 1, 1000]}
	I_units_df = pd.DataFrame(I_units_df)
	I_units_df.index = ['A', 'mA', 'uA']

	# Define matrix for unit conversion for voltage
	I_units_df = {'A': [1, 60, 3600, 3600000], 'mA': [1/60, 1, 60, 60000], 'uA': [1/3600, 1/60, 1, 1000]}
	I_units_df = pd.DataFrame(I_units_df)
	I_units_df.index = ['A', 'mA', 'uA']

	# Define matrix for unit conversion for capacity
	C_units_df = {'Ah': [1, 60, 3600, 3600000], 'mAh': [1/60, 1, 60, 60000], 'uAh': [1/3600, 1/60, 1, 1000]}
	C_units_df = pd.DataFrame(C_units_df)
	C_units_df.index = ['Ah', 'mAh', 'uAh']

	







