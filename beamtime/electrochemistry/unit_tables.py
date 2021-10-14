import pandas as pd

def time():
    # Define matrix for unit conversion for time
    time = {'h': [1, 60, 3600, 3600000], 'min': [1/60, 1, 60, 60000], 's': [1/3600, 1/60, 1, 1000], 'ms': [1/3600000, 1/60000, 1/1000, 1]}
    time = pd.DataFrame(time)
    time.index = ['h', 'min', 's', 'ms']

    return time

def current():
    # Define matrix for unit conversion for current
    current = {'A': [1, 1000, 1000000], 'mA': [1/1000, 1, 1000], 'uA': [1/1000000, 1/1000, 1]}
    current = pd.DataFrame(current)
    current.index = ['A', 'mA', 'uA']

    return current

def voltage():
    # Define matrix for unit conversion for voltage
    voltage = {'V': [1, 1000, 1000000], 'mV': [1/1000, 1, 1000], 'uV': [1/1000000, 1/1000, 1]}
    voltage = pd.DataFrame(voltage)
    voltage.index = ['V', 'mV', 'uV']

    return voltage

def capacity():
    # Define matrix for unit conversion for capacity
    capacity = {'Ah': [1, 1000, 1000000], 'mAh': [1/1000, 1, 1000], 'uAh': [1/1000000, 1/1000, 1]}
    capacity = pd.DataFrame(capacity)
    capacity.index = ['Ah', 'mAh', 'uAh']

    return capacity

def mass():
    # Define matrix for unit conversion for capacity
    mass = {'kg': [1, 1000, 1000000, 1000000000], 'g': [1/1000, 1, 1000, 1000000], 'mg': [1/1000000, 1/1000, 1, 1000], 'ug': [1/1000000000, 1/1000000, 1/1000, 1]}
    mass = pd.DataFrame(mass)
    mass.index = ['kg', 'g', 'mg', 'ug']

    return mass


def energy():
    
    energy = {'kWh': [1, 1000, 1000000], 'Wh': [1/1000, 1, 1000], 'mWh': [1/100000, 1/1000, 1]}
    energy = pd.DataFrame(energy)
    energy.index = ['kWh', 'Wh', 'mWh']

    return energy



