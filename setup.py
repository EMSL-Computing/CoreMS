import pathlib, os, sys

import setuptools
from setuptools import Command, setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

class test_enviroms(Command):
    
    """Run all of the tests for the package.
    This is a automatic test run class to make distutils 
    python setup.py build
    python setup.py install
    python setup.py test
    """

    description = "Automatically run the test suite for Biopython."
    user_options = [("offline", None, "Don't run online tests")]

    def initialize_options(self):
        """No-op, initialise options."""
        #self.offline = None
        pass

    def finalize_options(self):
        """No-op, finalise options."""
        pass

    def run(self):
        """Run the tests."""
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("tests")
        sys.path.insert(0, "")
        import run_test
        #if self.offline:
        #    run_test.main(["--offline"])
        #else:
        run_test.main([])

        # change back to the current directory
        os.chdir(this_dir)
        
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
    packages= setuptools.find_packages(".", exclude= ["test", "*win_only"]),
    exclude_package_data={'.': ["test", "*.win_only"]},
    include_package_data=True,
    install_requires=["pandas", "numpy", "matplotlib", "scipy", 'IsoSpecPy'],
    # test are not yet implemented, will test dependences and syntax only for now
    test_suite='pytest',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    #test_suite='nose2.collector',
    #tests_require=['nose2'],
    #cmdclass={
          #"install": install_enviroms,
          #"build_py": build_py_enviroms,
          #"build_ext": build_ext_enviroms,
    #      "test": test_enviroms,
    #  },
    #entry_points={
    #    "console_scripts": [
    #        "enviroms=cli.__main__:main",
    #    ]
    #},
)
