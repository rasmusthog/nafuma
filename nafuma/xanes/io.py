import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import nafuma.auxillary as aux 


def split_xanes_scan(root, destination=None, replace=False):
    #root is the path to the beamtime-folder
    #destination should be the path to the processed data
    
    #insert a for-loop to go through all the folders.dat-files in the folder root\xanes\raw

    # FIXME Only adding this variable to pass the Linting-tests - will refactor this later
    filename = 'dummy'
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    datas = []
    data = []
    headers = []
    header = ''
    start = False
    
    for line in lines:
        if line[0:2] == "#L":
            start = True
            header = line[2:].split()
            continue
            
        elif line[0:2] == "#C":
            start = False
            
            if data:
                datas.append(data)
                data = []
                
            if header:
                headers.append(header)
                header = ''
                
                

        if start == False:
            continue
            
        else:
            data.append(line.split())
            
            
            
            
    edges = {'Mn': [6.0, 6.1, 6.2, 6.3, 6.4, 6.5], 'Fe': [6.8, 6.9, 7.0, 7.1, 7.2], 'Co': [7.6, 7.7, 7.8, 7.9], 'Ni': [8.1, 8.2, 8.3, 8.4, 8.5]}
    edge_count = {'Mn': 0, 'Fe': 0, 'Co': 0, 'Ni': 0}
    

    for ind, data in enumerate(datas):
        df = pd.DataFrame(data)
        df.columns = headers[ind]

        edge_start = np.round((float(df["ZapEnergy"].min())), 1)

        for edge, energies in edges.items():
            if edge_start in energies:
                edge_actual = edge
                edge_count[edge] += 1

        
        
        filename = filename.split('/')[-1]
        count = str(edge_count[edge_actual]).zfill(4)

        
        # Save 
        if destination:
            cwd = os.getcwd()

            if not os.path.isdir(destination):
                os.mkdir(destination)
                
            os.chdir(destination)

            df.to_csv('{}_{}_{}.dat'.format(filename.split('.')[0], edge_actual, count))

            os.chdir(cwd)

        else:
            df.to_csv('{}_{}_{}.dat'.format(filename.split('.')[0], edge_actual, count))






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
