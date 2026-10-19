FROM stewartlab/scalinger:latest

WORKDIR /src

# **Copy the STAR.sh script into /src**
RUN rm -rf ./src/STAR.sh
COPY STAR.sh /src/

# **Copy the GEM-X (v4) barcode whitelists alongside the V2/V3 ones**
COPY 3M-3pgex-may-2023_TRU.txt 3M-5pgex-jan-2023.txt /usr/local/bin/

# **Ensure that STAR.sh is executable**
RUN chmod +x /src/STAR.sh
