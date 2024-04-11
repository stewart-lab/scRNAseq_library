# Start with a Python base image
FROM jfreeman88/scrna-seq:bioinfo_latest
USER ${username:-$(id -u):$(id -g)}
ARG username
#RUN useradd -ms /bin/bash $(id -u):$(id -g)
#
USER $username
# Set the working directory
WORKDIR /scRNA-seq
# Copy the current directory contents into the container at /scRNA-seq
COPY --chown=$username . /scRNA-seq

