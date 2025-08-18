FROM stewartlab/scrnaseq-env:latest
RUN chmod -R 777 /root
RUN apt-get update && apt-get install -y locales \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
RUN rm -rf /scRNA-seq
COPY ./ /scRNA-seq
WORKDIR /scRNA-seq