import pandas as pd
import numpy as np
import os

def rbkerbest():
    print("ROSENBORG!<3")


def split_xanes_scan(filename, destination=None):

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