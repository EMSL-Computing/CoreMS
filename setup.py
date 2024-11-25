import os
import pathlib
import sys
from setuptools import setup, find_packages
# The directory containing this file

# The text of the README file
with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

# This call to setup() does all the work
setup(
    name="CoreMS",
    version="3.1.0",
    description="Mass Spectrometry Framework for Small Molecules Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EMSL-Computing/CoreMS",
    author="Corilo, Yuri",
    author_email="corilo@pnnl.gov",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Development Status :: 4 - Beta"
    ],


    # package_data={'external': ['disclaimer.txt'], '': ['ext_lib/*']},
    packages=find_packages(),
    exclude_package_data={'.': ["tests", "*.win_only"]},
    include_package_data=True,
    install_requires=required,
    setup_requires=['pytest-runner', 'wheel'],
    test_suite='pytest',
    tests_require=['pytest'],
)
