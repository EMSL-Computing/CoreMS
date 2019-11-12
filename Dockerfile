FROM gitlab.pnnl.gov:4567/mass-spectrometry/corems:corems-base

WORKDIR /home/CoreMS

COPY corems/ /home/CoreMS/corems/
COPY requirements.txt LICENSE README.md setup.py /home/CoreMS/
RUN python3 -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('win_only/__init__.py')]"
RUN python3 -m pip install -r requirements.txt
RUN python3 setup.py install 
RUN pip3 install jupyter
CMD jupyter notebook --ip 0.0.0.0 --no-browser --allow-root


