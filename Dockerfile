FROM python:3.7.5-slim-buster
WORKDIR /home/corems/
#VOLUME  . /home/corems/

COPY corems/ /home/corems/
COPY requirements.txt LICENSE README.md setup.py /home/corems/
RUN apt install g++
RUN python3 -m venv corems_env
RUN source corems_env/bin/activate
RUN python3 -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('win_only/__init__.py')]"
RUN python3 -m pip install --upgrade setuptools wheel
RUN python3 -m pip install -r requirements.txt
RUN python3 -m pip install dist/*
#RUN deploy to PyPi 
#RUN build cli
#RUN deploy to dockerhub 


