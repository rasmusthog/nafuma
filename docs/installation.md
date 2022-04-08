# Installation

This package is not available on any package repositories, but can be installed by cloning the repository from GitHub and installing via ```pip install``` from the root folder:

```
$ git clone git@github.com:rasmusthog/nafuma.git
$ cd nafuma
$ pip install .
``` 
If you are planning on making changes to the code base, you might want to consider installing it in develop-mode in order for changes to take effect without reinstalling by including the ```-e``` flag: 

```
pip install -e .
```

As of now (v0.2, 08-04-22), the installer will not install any dependencies. It is recommended that you use `conda` to create an environment from `environment.yml` in the root folder:

```
$ conda env create --name <your_environment_name_here> --file environment.yml
$ conda activate <your_environment_name_here>
```

(remember to also get rid of <> when substituting your environment name).

This should get you up and running!
