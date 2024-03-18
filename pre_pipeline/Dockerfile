# Use the mambaforge base image
FROM jfreeman88/get_data:STAR

# Initialize shell for conda
RUN conda init bash

# Set working directory
WORKDIR /src

# Copy the current directory contents into the container at /src
COPY . /src/

# change permissions of everything in the /data directory
RUN chmod -R 777 /data/

RUN chmod +x /src/STAR.sh















