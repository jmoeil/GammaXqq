#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_103 x86_64-centos7-gcc12-opt
cd /user/jmoeil/GammaXqq/macros
python3 analysis.py --max_events -1 -i $1 -o $2 -p $3 --year $4 --era $5 --isData $6


