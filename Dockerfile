# Start with a Python base image
FROM jfreeman88/scrna-seq:latest

# Set the working directory
WORKDIR /scRNA-seq

COPY . /scRNA-seq

