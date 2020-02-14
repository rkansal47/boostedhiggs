#!/usr/bin/env bash

source env_lcg.sh

pip install --user coffea
pip install --upgrade coffea --user

# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension
