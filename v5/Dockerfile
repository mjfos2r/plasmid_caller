FROM python:3.9-slim

LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/plasmid_caller" \
    org.opencontainers.image.description="Python3 application for the classification of plasmids in Borrelia burgdorferi genome assemblies" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="mfoster11@mgh.harvard.edu"

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /tmp/blast /db /app
WORKDIR /app

RUN wget -P /tmp/blast https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
RUN wget -P /tmp/blast https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz.md5
RUN cd /tmp/blast && md5sum -c /tmp/blast/ncbi-blast-2.16.0+-x64-linux.tar.gz.md5 || exit 1

RUN tar -xzvf /tmp/blast/ncbi-blast-2.16.0+-x64-linux.tar.gz -C /tmp/blast/ --strip-components=1 --no-same-owner && \
    mv /tmp/blast/bin/* /usr/local/bin/
    
RUN rm -rf /tmp/blast

COPY ["db", "/db/"]
COPY ["requirements.txt", "plasmid_caller.py", "blast_parsing_dict.pkl", "/app/"]

RUN pip install --no-cache-dir -r /app/requirements.txt

ENTRYPOINT ["python", "/app/plasmid_caller.py"]
