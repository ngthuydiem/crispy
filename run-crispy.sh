#!/bin/zsh

#input: $INPUT, $K_CUTOFF, $CUTOFF, $B
INPUT=$1
K=$2
K_CUTOFF=$3
B=$4
CUTOFF=$5

UNIQUE_INPUT=$INPUT"_Clean"

./crispy-cuda/bin/preprocess -i $INPUT
./crispy-cuda/bin/kmerDist -i $UNIQUE_INPUT -k $K -t $K_CUTOFF
./crispy-cuda/bin/genDist -i $UNIQUE_INPUT -b $B -t $CUTOFF 

NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
NUM_FILES=`ls $UNIQUE_INPUT".ndist"* | wc -l`	
./crispy-cuda/bin/aveclust -i $UNIQUE_INPUT -n $NUM_READS -f $NUM_FILES -e $CUTOFF 

