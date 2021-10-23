# beamtime
A package for processing and analysis of data from beamtime at SNBL


# 1.0 Installation

### 1.1 Virtual environment

It is recommended to use a new dedicated virtual environment to use this package to ensure that all depedency versions are the same. This can be done by first creating a new environtment with the Python version 3.9.7 as follows:

```
conda create --name beamtime python=3.9.7
```

Here you can replace `beamtime` with any name you prefer. If you do, make sure you replace it in the subsequent commands as well. 

In order to use the virtual environment, you need to activate it. This is done by typing

```
conda activate beamtime
```

Note that you might have to initialise your shell first to allow for this. If so, you'll get a message saying so and you need to follow the instructions you get.

Once you are able to activate the virtual environment, you need to install all the packages required. These are listed in the requirements.txt file in the beamtime package folder, and can be installed by running the command (inside that folder):

```
conda install --file requirements.txt
```


Lastly, you might want to install this virtual environment as its own Jupyter kernel if you are planning to use Jupyter Notebook / Labs with this package. This way you won't have to activate the environment everytime you use it, you just create a Jupyter Notebook with this kernel. 

To do so, run the command:

```
python -m ipykernel install --user --name beamtime --display-name "Python 3.9.7 (Beamtime)"
```

Here you need to change `beamtime` if you named it something else, and the display name can be whatever you want it to be. Note that you need to have the package `ipykernel` installed for this, but it should be installed from running the install command above.

### 1.2 Installation of the package

In order to also use `beamtime` package itself, you need to install it. It is not uploaded to any package manager, but it can be installed from the main folder containing the `setup.py` file. Run the following command in this folder:

```
pip install .
```

# 2.0 The `electrochemistry` module

The `electrochemistry` module allows to plot galvanostatic cycling data from BioLogic, Neware and Batsmall. 

General use:

```py
import beamtime.electrochemistry as ec

path = 'path/to/data/file'
options = {
  'x_vals': 'specific_capacity',
  'y_vals': 'voltage',
  'active_material_weight': 4.3
  }
  
cycles, fig, ax = ec.plot.plot_gc(path=path, kind='neware', options=options)

```

Note that no options needs to be specified, all options will have default values that should make somewhat sense. A comprehensive list of options will be updated later.

The return values from the `plot_gc` function are: 
- `cycles`, a list of lists containing the charge and discharge cycles of each cycle in the form of a `pandas` DataFrame
- `fig`, the `matplotlib.pyplot' Figure object to allow for any modifications after initial plotting.
- `ax`, the 'matplotlib'pyplot' Axes object to allow for any modification after initial plotting.

If these are not required, you can simply assign the return values to an underscore instead. However, omitting assignment will print all the DataFrames. 
