# Start with a Python base image
FROM jfreeman88/scrna-seq:latest

# Update the system and install necessary tools
RUN apt-get update && apt-get install -y g++ make wget

# Download and install SRA Toolkit
RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -vxzf sratoolkit.tar.gz && \
    SRATOOLKIT_DIR=$(ls -d */ | grep sratoolkit) && \
    mv $SRATOOLKIT_DIR/bin/* /usr/local/bin/ && \
    rm -r sratoolkit.tar.gz $SRATOOLKIT_DIR

# Download and install STAR
RUN wget --output-document STAR.tar.gz https://github.com/alexdobin/STAR/archive/master.tar.gz && \
    tar -vxzf STAR.tar.gz && \
    STAR_DIR=$(ls -d */ | grep STAR) && \
    cd $STAR_DIR/source && \
    make && \
    mv STAR /usr/local/bin/ && \
    cd ../../ && \
    rm -r STAR.tar.gz $STAR_DIR

# Set the working directory
WORKDIR /scRNA-seq

# Copy the current directory contents into the container at /scRNA-seq
COPY . /scRNA-seq

