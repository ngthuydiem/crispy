#!/bin/zsh

#input: $INPUT, $K, $B
PROFILER=(/usr/bin/time -f "%e %M")
INPUT=$1
K=$2
B=$3

UNIQUE_INPUT=$INPUT"_Clean"

./crispy-cuda/bin/preprocess -i $INPUT
$PROFILER ./crispy-cuda/bin/kmerDist -i $UNIQUE_INPUT -k $K 
$PROFILER ./crispy-cuda/bin/genDist -i $UNIQUE_INPUT -b $B

NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
NUM_FILES=`ls $UNIQUE_INPUT".ndist"* | wc -l`	
#$PROFILER ./crispy-cuda/bin/aveclust -i $UNIQUE_INPUT -n $NUM_READS -f $NUM_FILES
$PROFILER ./crispy-cuda/bin/hcluster -i $UNIQUE_INPUT -n $NUM_READS -f $NUM_FILES
