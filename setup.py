#!/usr/bin/env python3
"""
M6AC - MODTRAN6-based Atmospheric Correction for AVIRIS-NG
"""
from setuptools import setup, find_packages
import os

# Read long description from README
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='m6ac',
    version='0.1.0',
    author='Judy Northrop',
    author_email='your.email@example.com',
    description='MODTRAN6-based Atmospheric Correction for AVIRIS-NG',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/M6AC',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Image Processing',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'matplotlib>=3.3.0',
        'spectral>=0.22.0',
        'scikit-learn>=0.24.0',
        'tqdm>=4.60.0',
    ],
    extras_require={
        'geo': ['gdal>=3.0.0'],
        'kriging': ['pykrige>=1.6.0'],
        'dev': ['pytest>=6.2.0', 'pytest-cov>=2.12.0'],
        'docs': ['sphinx>=4.0.0', 'sphinx-rtd-theme>=0.5.0'],
    },
    entry_points={
        'console_scripts': [
            'm6ac-correct=m6ac.cli:main_correct',
            'm6ac-lut=m6ac.cli:main_lut',
        ],
    },
    include_package_data=True,
    package_data={
        'm6ac': ['data/*.json', 'data/*.txt'],
    },
)
