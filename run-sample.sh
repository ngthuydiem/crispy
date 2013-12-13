#!/bin/zsh

cd ..
INPUT="crispy-cuda/data/V2Mice_sample1000"
K=6
K_CUTOFF=0.5
B=20
CUTOFF=0.3

./crispy-cuda/run-crispy.sh $INPUT $K $K_CUTOFF $B $CUTOFF

