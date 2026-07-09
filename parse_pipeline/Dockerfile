FROM stewartlab/parse-biosciences:latest

WORKDIR /src

# jq is required by parse.sh to read config.json but is not present in the base image
RUN apt-get update && apt-get install -y --no-install-recommends jq \
    && rm -rf /var/lib/apt/lists/*

# anndata==0.13.0 (pinned in the base image) fails to import on Python <3.15:
# it calls `from typing_extensions import sentinel`, but typing_extensions only
# exposes the newer `Sentinel` class. Fixed upstream in 0.13.1.
RUN /opt/miniconda3/bin/pip install --no-cache-dir 'anndata==0.13.1'

# splitpipe only requires pandas>=2.2, so the base image resolved pandas 3.0.3.
# Pandas 3.0 made Arrow-backed strict-dtype string columns the default, which
# breaks splitpipe's combine step (assigning a float aggregate into what it
# reads as a string column raises TypeError instead of silently upcasting).
# Pin to the 2.x line splitpipe was actually built against.
RUN /opt/miniconda3/bin/pip install --no-cache-dir 'pandas>=2.2,<3'

# **Copy the parse.sh script into /src**
#RUN rm -rf ./src/parse.sh
COPY parse.sh parse_genome_build.sh /src/

# **Ensure that parse.sh and parse_genome_build.sh are executable**
RUN chmod +x /src/parse.sh /src/parse_genome_build.sh