import pandas as pd

def time():
    # Define matrix for unit conversion for time
    time = {'h': [1, 60, 3600, 3600000], 'min': [1/60, 1, 60, 60000], 's': [1/3600, 1/60, 1, 1000], 'ms': [1/3600000, 1/60000, 1/1000, 1]}
    time = pd.DataFrame(time)
    time.index = ['h', 'min', 's', 'ms']

    return time


