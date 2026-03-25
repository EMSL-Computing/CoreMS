FROM python:3.13-slim AS base
WORKDIR /home/corems

# Install .NET 8 runtime via official install script (avoids APT keyring SHA1 issue)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates libssl3 libkrb5-3 zlib1g && \
    curl -sSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh && \
    chmod +x /tmp/dotnet-install.sh && \
    /tmp/dotnet-install.sh --runtime dotnet --channel 8.0 --install-dir /usr/local/dotnet && \
    rm /tmp/dotnet-install.sh && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV DOTNET_ROOT=/usr/local/dotnet
ENV PATH="${PATH}:/usr/local/dotnet"
ENV PYTHONNET_RUNTIME=coreclr
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1

# Install the corems package (pyproject.toml defines all dependencies)
# gcc is needed to compile ms-entropy's Cython extension; purged afterwards to keep the image lean
COPY pyproject.toml README.md disclaimer.txt ./
COPY corems/ ./corems/
RUN apt-get update && apt-get install -y --no-install-recommends gcc python3-dev && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install --no-cache-dir . && \
    apt-get purge -y gcc python3-dev && apt-get autoremove -y && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    rm -rf corems/


FROM base AS build
WORKDIR /home/corems

COPY examples/notebooks/*.ipynb README.md ./


RUN python3 -m pip install --no-cache-dir jupyter

ENV TINI_VERSION=v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]


