#!/bin/bash

########################################
############# CSCI 2951-O ##############
########################################
# Update this file with instructions on how to compile your code
if ! which uv ; then 
	echo "Installing uv"
	curl -LsSf https://astral.sh/uv/install.sh | sh
fi	
uv venv p4_venv --python 3.9
source p4_venv/bin/activate
uv pip install -r requirements.txt
cd /course/cs2951o/ilog/CPLEX_Studio2211/cplex/python/3.9/x86-64_linux/ && python3 setup.py install
