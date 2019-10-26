FROM python:3.7.5
WORKDIR /home/corems/
#VOLUME  . /home/corems/
#python:3.7.5-alpine
COPY corems/ /home/corems/
COPY requirements.txt LICENSE README.md setup.py /home/corems/
RUN python -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('win_only/__init__.py')]"
RUN python -m pip install --upgrade setuptools wheel
#RUN python3 -m pip install -r requirements.txt
RUN python setup.py sdist
RUN python -m pip install dist/*
#RUN deploy to PyPi 
#RUN build cli
#RUN deploy to dockerhub 


