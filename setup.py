import pathlib, os, sys

import setuptools
from setuptools import Command, setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="coreMS",
    version="3.0.0-alpha",
    description="Object Oriented Mass Spectrometry ToolBox for Small Molecules",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://gitlab.pnnl.gov/corilo/corems/",
    author="Corilo, Yuri",
    author_email="corilo@pnnl.gov",
    license="GNU Affero General Public License v3.0",
    
    classifiers=[
        "License :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    packages= setuptools.find_packages(".", exclude= ["test", "*win_only"]),
    exclude_package_data={'.': ["test", "*.win_only"]},
    include_package_data=True,
    install_requires=["pandas", "numpy", "matplotlib", "scipy", 'IsoSpecPy', 'sqlalchemy', 'openpyxl', 'pymongo', 'psycopg2-binary'],
    # test are not yet implemented, will test dependences and syntax only for now
    setup_requires=['pytest-runner', 'wheel'],
    test_suite='pytest',
    tests_require=['pytest'],
    #entry_points={
    #    "console_scripts": [
    #        "corems=cli.__main__:main",
    #    ]
    #},
)
