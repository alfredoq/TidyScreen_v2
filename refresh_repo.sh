#!/bin/bash

environment=$(echo $CONDA_DEFAULT_ENV)

rsync -rtuvzhl --delete ./tidyscreen /home/fredy/anaconda3/envs/$environment/lib/python3.12/site-packages

echo "Environment: $environment updated"

