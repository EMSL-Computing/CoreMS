FROM python:3.7.4
WORKDIR /home/corems/
#VOLUME  . /home/corems/

COPY  corems/ /home/corems/
COPY  requirements.txt LICENSE README.md setup.py /home/corems/
RUN python -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('win_only/__init__.py')]"
RUN pip3 install -r requirements.txt
RUN pip3 install wheel
RUN python setup.py sdist bdist_wheel
RUN pip install dist/*
#RUN deploy to PyPi 
#RUN build cli
#RUN deploy to dockerhub 


