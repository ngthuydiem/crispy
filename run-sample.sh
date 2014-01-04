#!/bin/zsh

cd ..
INPUT="crispy-cuda/data/V2Mice_sample1000"
K=6
B=20

./crispy-cuda/run-crispy.sh $INPUT $K $B

