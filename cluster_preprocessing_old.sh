#!/bin/bash
set -e

# Cluster version for quantifying several samples in paralel.
# Use: bash cluster/preprocessing.sh /path/to/config/file.config
 
# Paths
config_file=$1
source $config_file

# Run quantification
# Create_index
cat ${docDir}samples_id.txt | while read in
do
	#in="BM_Monocyte_rep1"
	qsub -P AG -A arubio -l h_vmem=4G -v config_file=$config_file,sample_name=$in -N "${in}_prepro" ${srcDir}01.preprocessing.sh
done


