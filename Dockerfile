FROM python:3.13-slim AS base
WORKDIR /home/corems

# Install .NET 8 runtime (required for pythonnet)
RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates && \
    wget -q https://packages.microsoft.com/config/debian/12/packages-microsoft-prod.deb \
         -O /tmp/microsoft-prod.deb && \
    dpkg -i /tmp/microsoft-prod.deb && rm /tmp/microsoft-prod.deb && \
    apt-get update && \
    apt-get install -y --no-install-recommends dotnet-runtime-8.0 && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV PYTHONNET_RUNTIME=coreclr

# Install Python dependencies as a separate layer for better cache reuse
COPY requirements.txt ./
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install --no-cache-dir -r requirements.txt

# Install the corems package
COPY pyproject.toml README.md disclaimer.txt ./
COPY corems/ ./corems/
RUN pip install --no-deps --no-cache-dir . && \
    rm -rf corems/


FROM base AS build
WORKDIR /home/corems

COPY examples/notebooks/*.ipynb README.md disclaimer.txt ./
COPY examples/scripts ./examples/

RUN python3 -m pip install --no-cache-dir jupyter

ENV TINI_VERSION=v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]


