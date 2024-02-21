FROM node:16-buster

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y less awscli jq \
#
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*


RUN mkdir -p /build

COPY . /build

WORKDIR /build

RUN npm install -g hirenj/node-reciperunner

RUN npm install

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=

CMD /build/bin/run_docker_conversion.sh