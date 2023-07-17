# Start with a Python base image
FROM python:3.11.3

# Set the working directory
WORKDIR /scRNA-seq

# Install necessary utilities
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    bash

# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p $HOME/miniconda && \
    rm ~/miniconda.sh && \
    ~/miniconda/bin/conda clean -afy && \
    ln -s ~/miniconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". ~/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

COPY . /scRNA-seq

# Switch shell back to default
SHELL ["/bin/bash", "-c"]

# Make the script executable
RUN chmod +x install_venv.sh

# Execute the script
RUN ./install_venv.sh scrnaseq

# Set the default command to start a Bash shell
CMD ["/bin/bash"]
