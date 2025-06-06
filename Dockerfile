FROM python:3.13-slim

LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/plasmid_caller" \
    org.opencontainers.image.description="Python3 application for the classification of plasmids in Borrelia burgdorferi genome assemblies" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="mfoster11@mgh.harvard.edu"

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    bash \
    && rm -rf /var/lib/apt/lists/*

# add UV
ADD https://astral.sh/uv/0.7.6/install.sh /uv-installer.sh
RUN sh /uv-installer.sh && rm /uv-installer.sh
ENV PATH="/root/.local/bin/:$PATH"
RUN mkdir -p /app/plasmid_caller /data

# add seqkit
RUN curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz >/tmp/seqkit_linux_amd64.tar.gz \
    && curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz.md5.txt >/tmp/seqkit_linux_amd64.tar.gz.md5.txt

RUN [ "$(md5sum '/tmp/seqkit_linux_amd64.tar.gz' | awk '{print $1}')" = "$(awk '{print $1}' /tmp/seqkit_linux_amd64.tar.gz.md5.txt)" ]\
    && { echo -e "\033[0;32mChecksum valid!\033[0m"; tar -xzf /tmp/seqkit_linux_amd64.tar.gz -C /usr/local/bin --no-same-owner; rm /tmp/seqkit_linux_amd64.tar.gz; } \
    || { echo -e "\033[0;31mChecksum invalid!\033[0m"; exit 1; }

WORKDIR /app/plasmid_caller
COPY . .
WORKDIR /app
RUN uv venv
ENV PATH="/app/.venv/bin:$PATH"
RUN uv pip install plasmid_caller/
RUN plasmid_caller --version
WORKDIR /data
ENTRYPOINT [ "/bin/bash" ]