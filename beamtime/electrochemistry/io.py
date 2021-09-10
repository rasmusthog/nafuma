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

	