#!/bin/zsh

#INPUT="data/V2Mice_sample30000"
#UNIQUE_INPUT=$INPUT"_Clean"

STEP_SIZE=1000
THRESHOLD=1.0

START=2
END=2

for (( i=$START; i<=$END; i++ )) 
do
(( NUM_READS = $i * $STEP_SIZE ))
NAME="microbial_"$NUM_READS

INPUT="/DATA/grinder_data"/$NAME"-reads.fa"
OUTPUT="/DATA/grinder_genmat"/$NAME
UNIQUE_INPUT=$INPUT

echo "#########################  CRISPY-CUDA (GPU) #########################"

START_TIME=$(date +%s.%N)
#./bin/preprocess -i $INPUT 
./bin/kmerDist -i $UNIQUE_INPUT -k 2 -t 1.0

NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
echo $NUM_READS
NUM_ENTRIES=`./bin/genDist -i $UNIQUE_INPUT -t $THRESHOLD -b 1`
NUM_FILES=`ls $UNIQUE_INPUT".ndist"* | wc -l`	

./bin/aveclust -i $UNIQUE_INPUT -m $NUM_ENTRIES -n $NUM_READS -f $NUM_FILES -o $OUTPUT -e $THRESHOLD -s 0.1

END_TIME=$(date +%s.%N)
DIFF=$(echo " $END_TIME - $START_TIME " | bc)

(( SPARSITY = $NUM_ENTRIES * 1.0 / ($NUM_READS * ($NUM_READS - 1) / 2) ))
echo "Sparse matrix contains "$NUM_ENTRIES" entries in "$NUM_FILES" files with sparsity "$SPARSITY
echo "It took $DIFF seconds"

rm -f $INPUT"."*

done


