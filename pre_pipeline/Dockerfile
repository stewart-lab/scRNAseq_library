FROM stewartlab/scalinger:latest

WORKDIR /src

# **Copy the STAR.sh script into /src**
RUN rm -rf ./src/STAR.sh
COPY STAR.sh /src/

# **Ensure that STAR.sh is executable**
RUN chmod +x /src/STAR.sh