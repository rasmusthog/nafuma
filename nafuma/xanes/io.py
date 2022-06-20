import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import nafuma.auxillary as aux
from nafuma.xanes.calib import find_element


def split_scan_data(data: dict, options={}):
    

    required_options = ['save', 'save_folder', 'replace', 'add_rois']

    default_options = {
        'save': False,
        'save_folder': '.',
        'replace': False,
        'add_rois': False
    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)
    #root is the path to the beamtime-folder
    #destination should be the path to the processed data
    
    #insert a for-loop to go through all the folders.dat-files in the folder root\xanes\raw

    # FIXME Only adding this variable to pass the Linting-tests - will refactor this later

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

    all_scans = []
    
    for filename in data['path']:

        with open(filename, 'r') as f:
            lines = f.readlines()
            
        scan_datas, scan_data = [], []
        headers, header = [], ''
        read_data = False
        
        for line in lines:
            # Header line starts with #L - reads headers, and toggles data read-in on
            if line[0:2] == "#L":
                header, read_data = line[2:].split(), True
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

            if not ('xmap_roi00' in headers[i]) and (not 'xmap_roi01' in headers[i]):
                continue 

            
            edge = find_element({'xanes_data_original': xanes_df})
            edges[edge].append(xanes_df)
            
        
        if options['add']:
   
            added_edges = {'Mn': [], 'Fe': [], 'Co': [], 'Ni': []}
            for edge, scans in edges.items():
                if scans:
                    xanes_df = scans[0]

                    for i, scan in enumerate(scans):
                        if i > 0:

                            if 'xmap_roi00' in xanes_df.columns:
                                xanes_df['xmap_roi00'] += scan['xmap_roi00']
                            if 'xmap_roi01' in xanes_df.columns:
                                xanes_df['xmap_roi01'] += scan['xmap_roi01']

                    added_edges[edge].append(xanes_df)

            edges = added_edges
            
        if options['save']:
            if not os.path.isdir(options['save_folder']):
                os.makedirs(options['save_folder'])


            filename = os.path.basename(filename).split('.')[0]

            for edge, scans in edges.items():
                for i, scan in enumerate(scans):
                    count = '' if options['add'] else '_'+str(i).zfill(4)
                    path = os.path.join(options['save_folder'], f'{filename}_{edge}{count}.dat')
                    scan.to_csv(path)

        all_scans.append(edges)


    return all_scans



def read_data(data: dict, options={}) -> pd.DataFrame:


    # FIXME Handle the case when dataseries are not the same size

    required_options = []
    default_options = {

    }

    options = aux.update_options(options=options, required_options=required_options, default_options=default_options)

    columns = ['ZapEnergy']

    # Initialise DataFrame with only ZapEnergy-column
    xanes_data = pd.read_csv(data['path'][0])[['ZapEnergy']]

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]

    for filename in data['path']:
        columns.append(filename)

        scan_data = pd.read_csv(filename)

        scan_data = scan_data[[determine_active_roi(scan_data)]]
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
    
    if (scan_data['xmap_roi00'].iloc[0:100].mean() < scan_data['xmap_roi00'].iloc[-100:].mean()) and (scan_data['xmap_roi01'].iloc[0:100].mean() < scan_data['xmap_roi01'].iloc[-100:].mean()):
        if (scan_data['xmap_roi00'].max()-scan_data['xmap_roi00'].min()) > (scan_data['xmap_roi01'].max() - scan_data['xmap_roi01'].min()):
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
