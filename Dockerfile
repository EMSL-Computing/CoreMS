FROM python:3.7.5-slim-buster
WORKDIR /home/corems/
#VOLUME  . /home/corems/
#python:3.7.5-alpine
COPY corems/ /home/corems/
COPY requirements.txt LICENSE README.md setup.py /home/corems/
RUN apt update
RUN apt install -y build-essential
RUN python -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('win_only/__init__.py')]"
#RUN python3 -m pip install -r requirements.txt
RUN python setup.py install 



