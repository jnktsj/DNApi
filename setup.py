#!/usr/bin/env python3

from setuptools import setup, find_packages
from dnapilib import __version__

setup(name='DNApi',
      version=__version__,
      description='De novo adapter prediction (iterative) algorithm for small \
      RNA sequencing data',
      author='Junko Tsuji',
      author_email='junko.tsuji@umassmed.edu',
      url='https://github.com/jnktsj/DNApi',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
      ],
      keywords='Adapter Prediction',
      packages=find_packages(),
      scripts=[
          'dnapi.py',
          'utils/qual_offset.py',
          'utils/qual_trim.py',
          'utils/to_fasta.py'
      ]
      )
