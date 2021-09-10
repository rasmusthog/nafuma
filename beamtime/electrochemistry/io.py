import pandas as pd


def read_battsmall(path):
	''' Reads BATTSMALL-data into a DataFrame.

	Input:
	path (required): string with path to datafile

	Output:
	df: pandas DataFrame containing the data as-is, but without additional NaN-columns.'''

	df = pd.read_csv(path, skiprows=2, sep='\t')
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

	df = unit_conversion(df=df, t=t, C=C, I=I, U=U)

	return df



def test_print():
	print('IT WORKS!')


def unit_conversion(df, t, C, I, U):

	print('YAHOO!')

	C, m = C.split('/')

	# Get the units used in the data set
	t_prev = df.columns[0].split()[-1].strip('[]')
	U_prev = df.columns[1].split()[-1].strip('[]')
	I_prev = df.columns[2].split()[-1].strip('[]')
	C_prev, m_prev = df.columns[4].split()[-1].strip('[]').split('/')

	print(t_prev)

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
	df["TT [{}]".format(t_prev)] = df["TT [{}]".format(t_prev)] * t_units_df[t_prev].loc[t]
	df["U [{}]".format(U_prev)] = df["U [{}]".format(U_prev)] * U_units_df[U_prev].loc[U]
	df["I [{}]".format(I_prev)] = df["I [{}]".format(I_prev)] * I_units_df[I_prev].loc[I]
	df["C [{}/{}]".format(C_prev, m_prev)] = df["C [{}/{}]".format(C_prev, m_prev)] * (C_units_df[C_prev].loc[C] / m_units_df[m_prev].loc[m])

	df.columns = ['TT [{}]'.format(t), 'U [{}]'.format(U), 'I [{}]'.format(I), 'Z1', 'C [{}/{}]'.format(C, m), 'Comment']



	return df

	






