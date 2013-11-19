#!/bin/zsh

INPUT="data/V2Mice_sample1000"
UNIQUE_INPUT=$INPUT"_Clean"

THRESHOLD=1.0

START=1
END=10

for (( i=$START; i<=$END; i++ )) 
do

echo "#########################  CRISPY-CUDA (GPU) #########################"

START_TIME=$(date +%s.%N)
./bin/preprocess -i $INPUT 
./bin/kmerDist -i $UNIQUE_INPUT -k 6 -t 1.0

NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
NUM_ENTRIES=`./bin/genDist -i $UNIQUE_INPUT -t $THRESHOLD -b $i`
NUM_FILES=`ls $UNIQUE_INPUT".ndist"* | wc -l`	

./bin/aveclust -i $UNIQUE_INPUT -m $NUM_ENTRIES -n $NUM_READS -f $NUM_FILES -e $THRESHOLD -s 0.1

END_TIME=$(date +%s.%N)
DIFF=$(echo " $END_TIME - $START_TIME " | bc)

(( SPARSITY = $NUM_ENTRIES * 1.0 / ($NUM_READS * ($NUM_READS - 1) / 2) ))
echo "Sparse matrix contains "$NUM_ENTRIES" entries in "$NUM_FILES" files with sparsity "$SPARSITY
echo "It took $DIFF seconds"

python2.7 ../parse-esprit.py -i $UNIQUE_INPUT".Cluster"
python2.7 ../compute-accuracy.py $UNIQUE_INPUT $UNIQUE_INPUT".Cluster.Formatted"

rm -f $UNIQUE_INPUT $UNIQUE_INPUT"."*

done


