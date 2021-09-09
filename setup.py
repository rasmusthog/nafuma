from setuptools import setup, find_packages

setup(name='beamtime',
      version='0.1',
      description='Package for process and analysis of beamtime data from SNBL',
      url='http://github.uio.no/rasmusvt/beamtime',
      author='Rasmus Vester Thøgersen, Halvor Høen Hval',
      author_email='rasmusvt@smn.uio.no',
      license='MIT',
      packages=find_packages(),
      zip_safe=False)
