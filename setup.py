import pathlib
from setuptools import setup
import setuptools


# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="enviroms",
    version="0.1.0",
    description="Object Oriented Mass Spectrometry ToolBox for small molecules",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://gitlab.pnnl.gov/corilo/enviroms/",
    author="Corilo, Yuri",
    author_email="corilo@pnnl.gov",
    license="Not decided yet",
    classifiers=[
        #"License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages= setuptools.find_packages(".", exclude= ["example", "*win_only"]),
    exclude_package_data={'.': ["example", "*.win_only"]},
    include_package_data=True,
    install_requires=["pandas", "numpy", "matplotlib", "scipy", 'IsoSpecPy'],
    #entry_points={
    #    "console_scripts": [
    #        "enviroms=cli.__main__:main",
    #    ]
    #},
)

