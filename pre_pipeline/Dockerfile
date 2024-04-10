# Use the mambaforge base image
FROM jfreeman88/get_data:STAR

ARG USER_ID
ARG GROUP_ID
RUN groupadd -g ${GROUP_ID} myuser && \
    useradd -l -u ${USER_ID} -g myuser -m myuser


# Initialize shell for conda
RUN conda init bash

# Set working directory
WORKDIR /src

# Copy the current directory contents into the container at /src
COPY . /src/

# change permissions of everything in the /data directory
RUN chmod -R 777 /data/

RUN chmod +x /src/STAR.sh















