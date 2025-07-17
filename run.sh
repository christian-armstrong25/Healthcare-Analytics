#!/bin/bash

########################################
############# CSCI 2951-O ##############
########################################
E_BADARGS=65
if [ $# -ne 1 ]
then
	echo "Usage: `basename $0` <input>"
	exit $E_BADARGS
fi

input=$1
source p4_venv/bin/activate || exit 1
# change this to point to your local installation
# CHANGE it back to this value before submitting
export DOCPLEX_COS_LOCATION=/course/cs2951o/ilog/CPLEX_Studio2211/CPLEX_Studio2211
# run the solver
python3.9 src/main.py $input
