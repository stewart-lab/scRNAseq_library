# Start with a Python base image
FROM stewartlab/seuratv5:1.01

# Create a user with a specific UID and GID, and create the home directory
ARG USER_ID
ARG GROUP_ID
RUN groupadd -g ${GROUP_ID} myuser && \
    useradd -l -u ${USER_ID} -g myuser -m myuser

# Set the working directory
WORKDIR /scRNA-seq

# Copy the current directory contents into the container at /scRNA-seq
COPY . /scRNA-seq

# Change file ownership to the new user
RUN chown -R myuser:myuser /scRNA-seq

RUN chown -R myuser:myuser /root
