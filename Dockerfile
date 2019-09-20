FROM python:3
WORKDIR /home/enviroms/
#VOLUME  . /home/enviroms/

#COPY  . /home/enviroms/
RUN   python -V  # Print out python version for debugging
RUN   pip install virtualenv
RUN   virtualenv venv
RUN   source venv/bin/activate
RUN   pip install -r requirements.txt
RUN   pip install wheel
RUN   python setup.py sdist bdist_wheel
RUN   python setup.py bdist_wheel
RUN   pip install dist/*