"""
"""

from setuptools import setup, find_packages
from codecs import open
from os import path



here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='rxnnet',
      version='0.1.1',
      description=('A Python package for representing reaction networks '
                   'and simulating their behaviors'),
      long_description=long_description,
      packages=find_packages(exclude=['tests', 'docs']),
      url='https://github.com/leihuang/rxnnet',
      license='MIT',
      author='Lei Huang',
      author_email='lh389@cornell.edu',
      description='kinetic modeling package',
      install_requires=['SloppyCell',
                        'pandas',
                        'numpy',
                        'matplotlib',
                        'scipy',
                        ],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Programming Language :: Python :: 2.7'],
)
