#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

podman build -t $CONF_project_docker - < $CONF_project_root_host/project.Dockerfile
