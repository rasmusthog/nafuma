import pandas as pd
import matplotlib.pyplot as plt
import os


def split_xanes_scan(root, destination=None, replace=False):
    #root is the path to the beamtime-folder
    #destination should be the path to the processed data
    
    #insert a for-loop to go through all the folders.dat-files in the folder root\xanes\raw
    
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


#Function that "collects" all the files in a folder, only accepting .dat-files from xanes-measurements
def get_filenames(path):
    
    
    cwd = os.getcwd()
    
    # Change into path provided
    os.chdir(path)
    
    filenames = [os.path.join(path, filename) for filename in os.listdir() if os.path.isfile(filename) and filename[-4:] == '.dat'] #changed
    
    
    
    # Change directory back to where you ran the script from
    os.chdir(cwd)
    
    return filenames

def put_in_dataframe(path):
    filenames = get_filenames(path) 

    #making the column names to be used in the dataframe, making sure the first column is the ZapEnergy
    column_names = ["ZapEnergy"]

    for i in range(len(filenames)):
        column_names.append(filenames[i])

    #Taking the first file in the folder and extracting ZapEnergies and intensity from that (only need the intensity from the rest)
    first = pd.read_csv(filenames[0], skiprows=0)

    #Making a data frame with the correct columns, and will fill inn data afterwards
    df = pd.DataFrame(columns = column_names)
    #First putting in the 2theta-values
    df["ZapEnergy"]=first["ZapEnergy"]

    #filling in the intensities from all files into the corresponding column in the dataframe
    for i in range(len(filenames)):
        df2 = pd.read_csv(filenames[i])
        df2 = df2.drop(['Mon','Det1','Det2','Det3','Det4','Det5', 'Det6','Ion1'], axis=1) #, axis=1)
        df2 = df2.drop(['MonEx','Ion2','Htime','MusstEnc1','MusstEnc3','MusstEnc4', 'TwoTheta', 'ZCryo'], axis=1)
        df2 = df2.drop(['ZBlower1', 'ZBlower2', 'ZSrcur'], axis=1)#, axis=19) #removing the sigma at this point
        
    ##############     THIS PART PICKS OUT WHICH ROI IS OF INTEREST, BUT MUST BE FIXED IF LOOKING AT THREE EDGES (roi00,roi01,roi02)    #####################
        if 'xmap_roi01' in df2.columns: 
            #Trying to pick the roi with the highest difference between maximum and minimum intensity --> biggest edge shift
            if max(df2["xmap_roi00"])-min(df2["xmap_roi00"])>max(df2["xmap_roi01"])-min(df2["xmap_roi01"]):
                df[filenames[i]]=df2["xmap_roi00"] #forMn
            else: 
                df[filenames[i]]=df2["xmap_roi01"] #forNi
        else:
            df[filenames[i]]=df2["xmap_roi00"]
    ###############################################################################################

        i=i+1


    #print(df)
    #If I want to make a csv-file of the raw data. Decided that was not necessary:
    #df.to_csv('static-Mn-edge.csv') #writing it to a csv, first row is datapoint (index), second column is 2theta, and from there the scans starts


    return df