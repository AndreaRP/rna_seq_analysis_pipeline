#!/bin/bash
set -e

# Cluster version for quantifying several samples in paralel.
# Use: bash cluster/quantification.sh /path/to/config/file.config
 
# Paths
config_file=$1
source $config_file

# Run quantification
# Create_index
cat ${docDir}samples_id.txt | while read in
do
	#in="BM_Monocyte_rep1"
	qsub -P AG -A arubio -l h_vmem=4G -v config_file=$config_file,index_file="mus_musculus_GRCm38.idx",sample_name=$in -N "${in}_quantification" ${srcDir}02.quantification.sh
done


