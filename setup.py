#!/usr/bin/env python
#  """HLA-Genotyper: Predict HLA genotypes from RNA-Seq and DNA-Seq Data."""
#
# from distutils.core import setup
from os import path
# from ez_setup import use_setuptools
# use_setuptools()
from setuptools import setup,find_packages

from codecs import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

execfile("hla/version.py")
setup(name='hla-genotyper',
      description='Tool to predict HLA genotypes from RNA-Seq and DNA-Seq Data',
      long_description=long_description,  
      version=version,
      author='John J. Farrell',
      author_email='farrell@bu.edu',
      url='http://github.com/jjfarrell/hla-genotyper',
      packages=find_packages(exclude=['docs','tests']),
      package_data={'hla': ['data/*.dat','data/*.txt',]},
      entry_points={
              'console_scripts': ['hla-genotyper = hla:main',],},
      classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",],
      install_requires=['numpy',
                        'pysam',
                        'biopython',
                       ],
      license='MIT',
      zip_safe=False)
