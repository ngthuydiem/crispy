#input: $INPUT, $K, $B
INPUT=$1
K=$2
B=$3

UNIQUE_INPUT=$INPUT"_Clean"

if [ -z "$INPUT" ]; then
	echo "Correct syntax: ./crispy.sh <input file name> [k value] [band value]"
else
	if [ -z "$K" ]; then 
		K=6
	fi
	if [ -z "$B" ]; then 
		B=10
	fi

	./bin/preprocess -i $INPUT
	./bin/kmerDist -i $UNIQUE_INPUT -k $K 
	./bin/genDist -i $UNIQUE_INPUT -b $B

	NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
	NUM_FILES=`ls $UNIQUE_INPUT".ndist"* | wc -l`	
	./bin/aveclust -i $UNIQUE_INPUT -n $NUM_READS -f $NUM_FILES
fi
