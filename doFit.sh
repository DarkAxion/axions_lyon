#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "Illegal number of parameters"
    echo "./doFit.sh [fitMin] [fitMax] [doBkg] [doAxions] [fixAxion]"
    exit 1
fi

fitMin=$1
fitMax=$2
doBkg=$3
doAxion=$4
fixAxion=$5

if [ $3 = 1 ]; then
  cd ../MC
  root make_histos.C+
  cd ../axions_lyon
fi


if [ $4 = 1 ]; then
  root make_model.C+
fi

####################### Perform the fit #############################

root fitter.C+() 
