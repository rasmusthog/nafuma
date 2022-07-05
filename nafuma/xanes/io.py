import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import nafuma.auxillary as aux
from nafuma.xanes.calib import find_element
from datetime import datetime


def split_scan_data(data: dict, options={}) -> list:
    ''' Splits a XANES-file from BM31 into different files depending on the edge. Has the option to add intensities of all scans of same edge into the same file. 
    As of now only picks out xmap_rois (fluoresence mode) and for Mn, Fe, Co and Ni K-edges.'''
    

    required_options = ['log', 'logfile', 'save', 'save_folder', 'replace', 'active_roi', 'add_rois', 'return']

    default_options = {
        'log': False,
        'logfile': f'{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}_split_edges.log',
        'save': False, # whether to save the files or not
        'save_folder': '.', # root folder of where to save the files
        'replace': False, # whether to replace the files if they already exist
        'active_roi': None,
        'add_rois': False, # Whether to add the rois of individual scans of the same edge together
        'return': True
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
            
        scan_datas, scan_data = [], []
        headers, header = [], ''
        read_data = False
        
        for line in lines:
            # Header line starts with #L - reads headers, and toggles data read-in on
            if line[0:2] == "#L":
                header, read_data = line[2:].split(), True

                if options['log']:
                    aux.write_log(message='... Found scan data. Starting read-in...', options=options)
                continue

            # First line after data started with #C - stops data read-in
            elif line[0:2] == "#C":
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
            
            xanes_df = pd.DataFrame(scan_data).apply(pd.to_numeric)
            xanes_df.columns = headers[i]
            edge = find_element({'xanes_data_original': xanes_df})

            if options['log']:
                aux.write_log(message=f'... Starting data clean-up ({edge}-edge)... ({i+1}/{len(scan_datas)})', options=options)


            if not ('xmap_roi00' in headers[i]) and (not 'xmap_roi01' in headers[i]):
                if options['log']:
                    aux.write_log(message='... ... Did not find fluoresence data. Skipping...', options=options)

                continue 

            
            
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
                        scan.to_csv(path)
                        if options['log']:
                            aux.write_log(message=f'... ... Scan saved to {path}', options=options)
                    
                    elif options['replace'] and os.path.isfile(path):
                        scan.to_csv(path)
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



def read_data(data: dict, options={}) -> pd.DataFrame:


    # FIXME Handle the case when dataseries are not the same size

    required_options = ['adjust']
    default_options = {
        'adjust': 0
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    columns = ['ZapEnergy']

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]
        
    # Initialise DataFrame with only ZapEnergy-column
    xanes_data = pd.read_csv(data['path'][0])[['ZapEnergy']]
    xanes_data['ZapEnergy'] += options['adjust']


    for filename in data['path']:
        columns.append(filename)

        scan_data = pd.read_csv(filename)

        if not options['active_roi']:
            scan_data = scan_data[[determine_active_roi(scan_data)]]
        else:
            scan_data = scan_data[options['active_roi']]
            
        xanes_data = pd.concat([xanes_data, scan_data], axis=1)


    xanes_data.columns = columns


    return xanes_data






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
