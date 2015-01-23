#!/bin/zsh
PROFILER=(/usr/bin/time -f "%e %M")
PROGRAM=(tools/crispy.sh)

# run the program
INPUT=$1
LOG=$2
SUMMARY=$3

if [[ -z $INPUT ]]; then
	echo 'Empty $INPUT!'
else
$PROFILER $PROGRAM $INPUT >>$LOG 2>>$SUMMARY

# postprocessing for accuracy evaluation
UNIQUE_INPUT=$INPUT"_Clean"
CLUSTER=$UNIQUE_INPUT"_Cluster"
FINAL_CLUSTER=$UNIQUE_INPUT"_FinalCluster"
MERGE=$UNIQUE_INPUT"_Merge"
FIG=$UNIQUE_INPUT"_Align.pdf"

R --no-save < crispy_changepoint.r --args $MERGE $FIG >>$LOG 2>>$LOG

NUM_READS=`grep '>' $UNIQUE_INPUT | wc -l`
./tools/crispy/bin/sparsecut -i $UNIQUE_INPUT -n $NUM_READS -c "optimal_cutoff.txt" >> $LOG

python parse_cluster.py -i $UNIQUE_INPUT".Cluster" -o $CLUSTER
#python compute_accuracy_crispy.py $UNIQUE_INPUT $CLUSTER >>$SUMMARY

python -W ignore find_clusters.py $UNIQUE_INPUT'.txt' $CLUSTER $FINAL_CLUSTER compute_accuracy_crispy.py $UNIQUE_INPUT $FINAL_CLUSTER >>$SUMMARY
python plots/tsne.py $UNIQUE_INPUT

# clean up	
rm 'optimal_cutoff.txt' $MERGE $FIG $CLUSTER #$FINAL_CLUSTER
rm $UNIQUE_INPUT $UNIQUE_INPUT'.npair_'* $UNIQUE_INPUT'.ndist_'*

EXTS=('.Cluster' '.Cluster_List' '.frq' '.kpair')
for EXT in ${EXTS[@]}; do
	FILE=$UNIQUE_INPUT$EXT
	rm $FILE
done


fi
