import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import nafuma.auxillary as aux
from nafuma.xanes.calib import find_element
import datetime

def split_scan_data(data: dict, options={}) -> list:
    ''' Splits a XANES-file from BM31 into different files depending on the edge. Has the option to add intensities of all scans of same edge into the same file. 
    As of now only picks out xmap_rois (fluoresence mode) and for Mn, Fe, Co and Ni K-edges.'''

    required_options = ['log', 'logfile', 'save', 'save_folder', 'replace', 'active_roi', 'add_rois', 'return', 'skip_if_no_roi']

    default_options = {
        'log': False,
        'logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_split_edges.log',
        'save': False, # whether to save the files or not
        'save_folder': '.', # root folder of where to save the files
        'replace': False, # whether to replace the files if they already exist
        'active_roi': None,
        'add_rois': False, # Whether to add the rois of individual scans of the same edge together
        'return': True,
        'skip_if_no_roi': True
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

    all_scans = []

    if options['log']:
        aux.write_log(message='Starting file splitting...', options=options)
    
    for filename in data['path']:

        if options['log']:
            aux.write_log(message=f'Reading {filename}...', options=options)

        with open(filename, 'r') as f:
            lines = f.readlines()
            
        timestamps = []
        scan_datas, scan_data = [], []
        headers, header = [], ''
        read_data = False
        
        for i, line in enumerate(lines):
            # Header line starts with #L - reads headers, and toggles data read-in on
            if 'zapline mono' in line:
                timestamps.append(lines[i+1].strip('#D'))
            
            elif line[0:2] == "#L":
                header, read_data = line[2:].split(), True

                if options['log']:
                    aux.write_log(message='... Found scan data. Starting read-in...', options=options)
                continue

            # First line after data started with #C - stops data read-in
            elif line[0:2] == "#C" or line[0:2] == '#S':
                read_data = False
                
                if scan_data:
                    scan_datas.append(scan_data); scan_data = []
                    
                if header:
                    headers.append(header); header = ''
                    
            # Ignore line if read-in not toggled       
            if read_data == False:
                continue
            
            # Read in data if it is
            else:
                scan_data.append(line.split())
                
                
        edges = {'Mn': [], 'Fe': [], 'Co': [], 'Ni': []}
        

        for i, scan_data in enumerate(scan_datas):

            if 'ZapEnergy' not in headers[i]:
                if options['log']:
                    aux.write_log(message=f'... No valid scan data found... ({i+1}/{len(scan_datas)})', options=options)
                continue
                        
            xanes_df = pd.DataFrame(scan_data).apply(pd.to_numeric)
            xanes_df.columns = headers[i]

            
            edge = find_element({'xanes_data_original': xanes_df})
        
                

            if options['log']:
                aux.write_log(message=f'... Starting data clean-up ({edge}-edge)... ({i+1}/{len(scan_datas)})', options=options)


            if not ('xmap_roi00' in headers[i]) and (not 'xmap_roi01' in headers[i]):
                if options['skip_if_no_roi']:
                    if options['log']:
                        aux.write_log(message='... ... Did not find fluoresence data. Skipping...', options=options)
                    continue 
                if options['log']:
                    aux.write_log(message='... ... Did not find fluoresence data, but still proceeding ...', options=options)

                    

            
            
            edges[edge].append(xanes_df)
            
        
        if options['add_rois']:

            if options['log']:
                aux.write_log(message=f'... Addition of rois enabled. Starting addition...', options=options)
   
            added_edges = {'Mn': [], 'Fe': [], 'Co': [], 'Ni': []}
            for edge, scans in edges.items():

                if options['log']:
                    aux.write_log(message=f'... ... Adding rois of the {edge}-edge...', options=options)
                
                if scans:
                    xanes_df = scans[0]

                    for i, scan in enumerate(scans):
                        if i > 0:

                            if options['log']:
                                aux.write_log(message=f'... ... ... Adding {i+1}/{len(scans)}', options=options)

                            if 'xmap_roi00' in xanes_df.columns:
                                xanes_df['xmap_roi00'] += scan['xmap_roi00']
                            if 'xmap_roi01' in xanes_df.columns:
                                xanes_df['xmap_roi01'] += scan['xmap_roi01']

                    added_edges[edge].append(xanes_df)

            edges = added_edges
            
        if options['save']:
            #FIXME If there is something wrong with the input file, the file will not be saved but log-file still sais it is saved. Goes from "Saving data to ..." to "All done!" no matter if it fals or not.
            if options['log']:
                aux.write_log(message=f'... Saving data to {options["save_folder"]}', options=options)

            if not os.path.isdir(options['save_folder']):
                if options['log']:
                    aux.write_log(message=f'... ... {options["save_folder"]} does not exist. Creating folder.', options=options)

                os.makedirs(options['save_folder'])


            filename = os.path.basename(filename).split('.')[0]

            for edge, scans in edges.items():
                for i, scan in enumerate(scans):
                    count = '' if options['add_rois'] else '_'+str(i).zfill(4)
                    path = os.path.join(options['save_folder'], f'{filename}_{edge}{count}.dat')
                    
                    if not os.path.isfile(path):
                        
                        with open(path, 'w', newline = '\n') as f:

                            f.write(f'# Time: {timestamps[i]}')
                            scan.to_csv(f)

                        if options['log']:
                            aux.write_log(message=f'... ... Scan saved to {path}', options=options)
                    
                    elif options['replace'] and os.path.isfile(path):
                        with open(path, 'w', newline = '\n') as f:
                            scan.to_csv(f)
                        
                        if options['log']:
                            aux.write_log(message=f'... ... File already exists. Overwriting to {path}', options=options)

                    elif not options['replace'] and os.path.isfile(path):
                        if options['log']:
                            aux.write_log(message=f'... ... File already exists. Skipping...', options=options)

        all_scans.append(edges)

    if options['log']:
        aux.write_log(message=f'All done!', options=options)


    if options['return']:
        return all_scans
    else:
        return



def save_data(data: dict, options={}) -> None:

    required_options = ['save_folder', 'overwrite', 'log', 'logfile', 'filename']

    default_options = {
        'log': False,
        'logfile': f'{datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_save_files.log',
        'save_folder': 'saved_scans',
        'overwrite': False,
        'filename': f'{datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_exported_data.dat',
        }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    # Check if there is any data to be saved
    if not 'xanes_data' in data.keys():
        if options['log']:
             aux.write_log(message=f'There is not saved scan data in data. Exiting without saving...', options=options)

        return None
    
    if not isinstance(data['xanes_data'], pd.DataFrame):
        if options['log']:
             aux.write_log(message=f'data["xanes_data"] has an invalid format. Exiting without saving...', options=options)

        return None


    # Make folder(s) if it/they do(es)n't exist
    if not os.path.exists(options['save_folder']):
        if options['log']:
             aux.write_log(message=f'Destination folder does not exist. Creating folder...', options=options)

        os.makedirs(options['save_folder'])



    if os.path.exists(os.path.join('save_folder', options['filename'])):
        if not options['overwrite']:
            if options['log']:
                aux.write_log(message=f'File already exists and overwrite disabled. Exiting without saving...', options=options)
            return None
        
    with open(os.path.join(options['save_folder'], options['filename']), 'w') as f:

        if 'e0_diff' in data.keys():
            f.write(f'# Number of header lines: {len(data["path"])+1} \n')

            for i, (path, e0) in enumerate(data['e0_diff'].items()):
                f.write(f'# Scan_{i} \t {e0} \n')
       
        else:
            f.write(f'# Number of header lines: {1}')
        

        data['xanes_data'].to_csv(f, sep='\t', index=False)


    #data['xanes_data'].to_csv(os.path.join(options['save_folder'], options['filename']), sep='\t', index=False)



def load_data(path: str) -> dict:
    # FIXME Let this function be called by read_data() if some criterium is passed

    data = {}


    with open(path, 'r') as f:
        line = f.readline()
        header_lines = int(line.split()[-1])

        if header_lines > 1:
            edge_positions = []
            line = f.readline()
            while line[0] == '#':
                edge_positions.append(line.split()[-1])
                line = f.readline()

    data['xanes_data'] = pd.read_csv(path, sep='\t', skiprows=header_lines)
    data['path'] = data['xanes_data'].columns.to_list()
    data['path'].remove('ZapEnergy')

    if header_lines > 1:
        data['e0_diff'] = {}

        for path, edge_position in zip(data['path'], edge_positions):
            data['e0_diff'][path] = float(edge_position)

    

    return data


def read_data(data: dict, options={}) -> pd.DataFrame:


    # FIXME Handle the case when dataseries are not the same size
    # FIXME Add possibility to extract TIME (for operando runs) and Blower Temp (for variable temperature runs)
    # FIXME Add possibility to iport transmission data
    required_options = ['adjust', 'mode']
    default_options = {
        'adjust': 0,
        'mode': 'fluoresence'
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    columns = ['ZapEnergy']

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]
        
    # Initialise DataFrame with only ZapEnergy-column
    xanes_data = pd.read_csv(data['path'][0], skiprows=1)[['ZapEnergy']]
    xanes_data['ZapEnergy'] += options['adjust']


    for filename in data['path']:
        columns.append(filename)

        scan_data = pd.read_csv(filename, skiprows=1)

        if options['mode'] == 'fluoresence':
            if not options['active_roi']:
                scan_data = scan_data[[determine_active_roi(scan_data)]]
            else:
                scan_data = scan_data[options['active_roi']]

        elif options['mode'] == 'transmission':
            scan_data = scan_data['MonEx'] / scan_data['Ion2']

        xanes_data = pd.concat([xanes_data, scan_data], axis=1)


    xanes_data.columns = columns


    return xanes_data


def read_metadata(data: dict, options={}) -> dict:

    required_options = ['get_temperature', 'get_timestamp', 'adjust_time', 'convert_time', 'time_unit', 'reference_time']

    default_options = {
        'get_temperature': True,
        'get_timestamp': True,
        'adjust_time': False,
        'convert_time': False,
        'reference_time': None,
        'time_unit': 's'
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)


    temperatures = []
    timestamps = []

    for filename in data['path']:
        scan_data = pd.read_csv(filename, skiprows=1)

        if options['get_temperature']:
            temperatures.append(scan_data['ZBlower2'].mean())

        if options['get_timestamp']:

            with open(filename, 'r') as f:
                #time = f.readline().strip('# Time: ') #<-- Previous code
                time = f.readline().split('# Time:  ')[-1] #Hope this does not fuck you up, Rasmus - but I needed another space here
                split_operator=time[-9] #This should be the operator that splits hours, minutes and seconds
                if split_operator == ".":
                    time = datetime.datetime.strptime(time, "%a %b %d %H.%M.%S %Y ")
                if split_operator == ":":
                    time = datetime.datetime.strptime(time, "%a %b %d %H:%M:%S %Y ")

            if options['adjust_time']:
                time_elapsed = scan_data['Htime'].iloc[-1] - scan_data['Htime'].iloc[0]

                time += datetime.timedelta(microseconds=time_elapsed)/2


            timestamps.append(time)


    if options['reference_time'] and options['convert_time']:
        from . import unit_tables
        new_times = []

        if isinstance(options['reference_time'], str):
            options['reference_time'] = datetime.datetime.strptime(options['reference_time'], "%d.%b %y %H.%M.%S")
        
        for time in timestamps:
            new_time = (time.timestamp() - options['reference_time'].timestamp()) * unit_tables.time()['s'].loc[options['time_unit']]        

            new_times.append(new_time)


        timestamps = new_times




    metadata = {'time': timestamps, 'temperature': temperatures}

    return metadata




def determine_active_roi(scan_data):

    # FIXME For Co-edge, this gave a wrong scan
    
    #Trying to pick the roi with the highest difference between maximum and minimum intensity --> biggest edge shift
    # if max(scan_data["xmap_roi00"])-min(scan_data["xmap_roi00"])>max(scan_data["xmap_roi01"])-min(scan_data["xmap_roi01"]):
    #     active_roi = 'xmap_roi00'
    # else: 
    #     active_roi = 'xmap_roi01'


    if not ('xmap_roi00' in scan_data.columns) or not ('xmap_roi01' in scan_data.columns):
        if 'xmap_roi00' in scan_data.columns:
            active_roi = 'xmap_roi00'
        elif 'xmap_roi01' in scan_data.columns:
            active_roi = 'xmap_roi01'
    
    elif (scan_data['xmap_roi00'].iloc[0:100].mean() < scan_data['xmap_roi00'].iloc[-100:].mean()) and (scan_data['xmap_roi01'].iloc[0:100].mean() < scan_data['xmap_roi01'].iloc[-100:].mean()):
        if (scan_data['xmap_roi00'].iloc[:int(scan_data.shape[0]/2)].max() - scan_data['xmap_roi00'].iloc[0])/scan_data['xmap_roi00'].max() > (scan_data['xmap_roi01'].iloc[:int(scan_data.shape[0]/2)].max() - scan_data['xmap_roi01'].iloc[0])/scan_data['xmap_roi01'].max():
            active_roi = 'xmap_roi00'
        else:
            active_roi = 'xmap_roi01'

    elif scan_data['xmap_roi00'].iloc[0:100].mean() < scan_data['xmap_roi00'].iloc[-100:].mean():
        active_roi = 'xmap_roi00'
    
    elif scan_data['xmap_roi01'].iloc[0:100].mean() < scan_data['xmap_roi01'].iloc[-100:].mean(): 
        active_roi = 'xmap_roi01'

    else:
        active_roi = None

    return active_roi
    


def write_data(data: dict, options={}):


    default_options = {
        'save_filenames': None,
        'save_dir': '.',
    }

    options = aux.update_options(options=options, default_options=default_options, required_options=default_options.keys())


    if not options['save_filenames']:
        options['save_filenames'] = [os.path.basename(col).split('.')[0]+'_exported.dat' for col in data['xanes_data'].columns if 'ZapEnergy' not in col]


    print(options['save_filenames'])