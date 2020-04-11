FROM corilo/corems:base-mono-pythonnet AS base
WORKDIR /home/CoreMS

COPY doc/notebooks/*.ipynb SettingsCoreMS.json /home/CoreMS/
COPY doc/examples /home/CoreMS/
COPY ESI_NEG_SRFA.d/ /home/CoreMS/ESI_NEG_SRFA.d/

#RUN apt update && apt install -y --no-install-recommends  build-essential

FROM base AS build  
COPY --from=base /home/CoreMS /home/CoreMS
WORKDIR /home/CoreMS

RUN python3 -m pip install corems
RUN python3 -m pip install jupyter
CMD jupyter notebook --ip 0.0.0.0 --no-browser --allow-root

