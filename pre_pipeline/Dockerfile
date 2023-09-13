# Use the mambaforge base image
FROM jfreeman88/get_data:STAR

# Initialize shell for conda
RUN conda init bash

# Set working directory
WORKDIR /src

# Copy the current directory contents into the container at /src
COPY . /src/

RUN chmod +x /src/STAR.sh















