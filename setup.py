"""
"""

from setuptools import setup, find_packages



setup(name='rxnnet',
      version='0.1.1',
      description=('A Python package for representing reaction networks '
                   'and simulating their behaviors'),
      packages=find_packages(exclude=['tests', 'docs']),
      url='https://github.com/leihuang/rxnnet',
      license='MIT',
      author='Lei Huang',
      author_email='lh389@cornell.edu',
      install_requires=['SloppyCell',
                        'pandas',
                        'numpy',
                        'matplotlib',
                        'scipy'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Programming Language :: Python :: 2.7'],
)
