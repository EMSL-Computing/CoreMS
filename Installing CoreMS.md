# Installing CoreMS
October 2023, Will Kew_

These instructions detail several different ways to get CoreMS up and running on your computer. Some scenarios may be simpler, while others may provide more development-functionality.  

As of CoreMS 2.0, minimum Python version supported is 3.10.


# Table of Contents
1. [CoreMS with Docker Database - Recommended](#medium)
2. [CoreMS with SQLite Database - Basic](#basic) 
3. [Additional Useful Installations](#misc)
4. [Development Installation - Advanced ](#dev)


## CoreMS with Docker Database - Recommended <a name="medium"></a>  
This approach is slightly more complex and recommended for those who want improved performance and the most up-to-date public releases of CoreMS. Here, we will be installing CoreMS from github, and using Docker to deploy a postgresql database.   

### Windows   
1. Download and install Python 3.10 
	a. [Python 3.10.11]](https://www.python.org/downloads/release/python-31011/)  
	b. *Windows installer (64-bit)*  
	c. During installation, make sure to select 'Add Python 3.10 to PATH'. **Note:** if you have another version of python installed, this may create a conflict.  
    d. Note - Python 3.10.11 is the last version of Py3.10 with a pre-built binary for Windows.   
2. Install Docker Desktop   
    a. https://www.docker.com/products/docker-desktop/  
3. [Optional] Install Git.   
	a.	 https://git-scm.com/downloads  
	b. If you do, you will more easily be able to update CoreMS in the future. If you dont, you will need to manually re-download CoreMS to get updates.   
4. [Git] Download CoreMS from Github   
	a. In a PowerShell prompt, navigate to a directory to download CoreMS to. For me, I create a 'Coding' folder for this.   
	b. `cd C:\`   
	c. `mkdir Coding`   
	d. `cd Coding`  
	e. `git clone https://github.com/EMSL-Computing/CoreMS.git`  
5. [Non-Git] Download CoreMS from Github  
	a. Navigate to https://github.com/EMSL-Computing/CoreMS  
	b. Click on 'Code' -> 'Download Zip'   
	c. Extract this to somewhere safe on your computer, e.g. C:\Coding
6. Set up a Python Virtual Environment to isolate CoreMS in.   
   From within a stable directory, i.e. within C:\Coding, execute:  
   a. `python -m venv venv`  
   b. `venv\Scripts\activate`    
   If this errors, you may need to allow scripts to be executed.  
   Launch PowerShell as an administrator and execute:  
   `set-executionpolicy RemoteSigned`    
   Then attempt to activate the venv again. 
   c. Note you can deactivate when finished with the similar command, `deactivate`  
   d. `pip install -U pip`     
   e. `pip install wheel setuptools`     
7. Install CoreMS with pip   
	a. Navigate to the directory you have put CoreMS, e.g. C:\Coding\corems, with a command prompt/powershell   
	b. `pip install .`   
	c. This should install CoreMS and any pre-requisites.   
8. For Thermo .raw file support, also install pythonnet  
    a. `pip install pythonnet`  
9.  Set up the PostgreSQL database using docker.   
    a. Within the corems directory, execute the command:  
    b. `docker-compose up -d`  
 
CoreMS should now be installed, and the docker database initialized. Note that you will need to be running Docker desktop to use this database, and it may not launch automatically on booting Windows.   

## MacOS
The MacOS Instructions may vary between versions of MacOS used, as well as specific needs of the user.   
Some of these steps (e.g. installed xcode-select and brew installing gcc, may take a significant amount of time.)  
1. Install command line tools:  
    a. `xcode-select --install`  
2. Install Homebrew on Mac  
    a. https://brew.sh/     
	b. Follow instructions on their website    
3. Install Python and Docker and Git and GCCwith Brew;  
   a. `brew install python@3.9`  
   b. `brew install git`  
   c. `brew install --cask docker`  
   d. `brew install gcc`  
4. Check the version of gcc installed and add a symlink for it    
   a. `find /usr/local/bin -name gcc*`      
   b. `ln -s g++-{version} g++`     
   c. `ln -s g++-{version} gcc`  
5.  If you are using Z shell, make sure to add the python path to your zprofile:    
   a. `echo 'eval "$(pyenv init --path)"' >> ~/.zprofile`    
6. Install Mono and DotNet. Necessary for handling Thermo Raw Files.     
   a. `brew install pkg-config`    
   b. `brew install mono-mdk`     
   c. `brew install mono`     
   d. `brew install dotnet`     
7. Clone the CoreMS repo to a location on your system which is 'stable', e.g.     
   a. `cd ~`     
   b. `mkdir CoreMS`    
   c. `cd CoreMS `    
   d. `git clone https://github.com/EMSL-Computing/CoreMS`    
8.  Set up a python virtual environment for CoreMS.    
    A virtual environment allows you to isolate the CoreMS packages and requirements from the rest of your python installation, providing more stability and robustness. It is strongly recommended that you use it.     
   a. `python3 -m venv venv`    
   b. `source venv/bin/activate`    
   c. `pip install -U pip`     
   d. `pip install wheel setuptools`    
9.  Find the installed version of Mono:  
   a. `find /Library/Frameworks/ -name mono-2.pc`  
   b. Run the following command with the found version of mono, e.g. 6.12.0  
   c. `$PKG_CONFIG_PATH=/Library/Frameworks/Mono.framework/Versions/{mono version}/lib/pkgconfig/ python3 -m pip install pythonnet==2.5.2`
10. Install CoreMS, CoreMS Dependencies  
    a. Enter the cloned CoreMS repo, e.g.   
	b. `cd corems`  
    c. `pip install -r requirements.txt`
11. Install Jupyter for browser based notebooks:  
    a. `pip install jupyter`
12. Install Spyder IDE:  
    a. `pip install spyder`
13. Spin up the Docker Database  
    a. `docker-compose up -d`


## CoreMS with SQLite Database - Basic <a name="basic"></a>  
This approach is the simplest to set up and will get you up and running fastest, but may not provide the best performance or most up-to-date features.  
This is based on the assumption you do not have Python 3 already installed. 
The major disadvantage to this method is that the molecular formula database will use sqlite and be slightly slower than the postgresql option, as well as being a local file rather than a system-wide accessible database.   
### Windows   
1. Download and install Python 3.10 
	a. [Python 3.10.11]](https://www.python.org/downloads/release/python-31011/) 
	b. *Windows installer (64-bit)*  
	c. During installation, make sure to select 'Add Python 3.9 to PATH'. **Note:** if you have another version of python installed, this may create a conflict.   
2. Launch a command prompt/powershell window and update and install key Python packages:  
	a.  `pip install --upgrade pip`   
	b.  `pip install --upgrade setuptools`   
	c.  `pip install wheel`   
3. Install CoreMS from PyPi, from within the command prompt.  
	a. `pip install corems`  
	b. If you want to be able to process Thermo .raw files, you'll also need to install PythonNet  
	c. `pip install pythonnet`  

The above steps should take less than 10 minutes and now you have CoreMS installed. Test the installation worked by loading up a Python window and importing the module.   
`import corems`  
`print(corems.__version__)`  
The output should be the version of CoreMS which was installed, e.g. 1.5.1  
	


## Additional Useful Installations <a name="misc"></a>  
In most cases, you will want to install or use an IDE to develop your CoreMS scripts and routines. Some options include Spyder, PyCharm, IDLE...   
I like to use Spyder for writing scripts and exploring data thanks its built in variable explorer.   
### Install Spyder with Pip  
`pip install spyder`  
This will get and install spyder. It has a number of dependencies so it may take a few minutes.   
You could also download Spyder with a standalone installer.  

### Install Jupyter Notebooks 
`pip install jupyter`  
This will download and install the requistes for Jupyter and Notebooks, including iPython if not already installed.  
You can then launch it from the command line with:
`jupyter notebook`  

### Install VisualStudio Code   
As CoreMS is a large and complex framework, it can be helpful to install a program for browsing the whole code base even as a non-developer. One good tool for this is VisualStudio Code.   
https://code.visualstudio.com/   
Download and install this, then 'Open Folder' and navigate to your CoreMS directory.   
If you installed CoreMS from PyPi, you will need to point it at the [python]/lib/site-packages/corems directory.   


## Development Installation - Advanced <a name="dev"></a>  
As above, but with additional installation of development packages for build testing.   

