#!/bin/bash
set -e

# Cluster version for quantifying several samples in paralel.
# Use: bash cluster/preprocessing.sh /path/to/config/file.config
while getopts "c:s:h" opt; do
  case $opt in
    h)
      echo "This will run the preprocessing script on the samples in DOC/samples_id.txt.
OPTIONS:
  -c: path to config file.
  -s: Single End (se) or Paired End (pe)
  Usage: bash ~/cluster_preprocessing.sh -c /path/to/config/file/file.config -s pe|se" 
      ;;
    c)
		config_file=$OPTARG
      ;;
    s)
    sequencing=$OPTARG
      ;;
    :)
    	echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done
# Paths
source $config_file

case $sequencing in
  pe)
    preprocessing_script="${srcDir}01.preprocessing_PE.sh"
    ;;
  se)
    preprocessing_script="${srcDir}01.preprocessing_SE.sh"
    ;;
  :)
    echo "Option $OPTARG requires an argument." 
    exit 1
    ;;
esac

array_length=($(wc -l ${docDir}samples_id.txt))

# Run quantification
cat ${docDir}samples_id.txt | while read in
do
	qsub -P AG -A arubio -l h_vmem=4G -v config_file=$config_file,sample_name=$in -N "${in}_prepro" $preprocessing_script
  # qsub -P AG -A arubio -l h_vmem=4G -N "${in}_prepro" -t 1-$array_length:2 ${srcDir}01.preprocessing.sh -c $config_file -s $in
done
