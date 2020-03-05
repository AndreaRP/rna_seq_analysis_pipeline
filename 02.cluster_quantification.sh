#!/bin/bash
set -e

# Cluster version for quantifying several samples in paralel.
# Use: bash cluster/preprocessing.sh /path/to/config/file.config
while getopts "c:q:t:h" opt; do
  case $opt in
    h)
      echo "This will run quantification on the processed reads.
  OPTIONS:
	-c: path to config file.
	-q: either 'rsem' or 'kallisto'
  -t: either single end (se) or paired end (pe)
	Usage: bash ~/cluster_quantification.sh -c /path/to/config/file/file.config -q rsem|kallisto -t se|pe" 
      ;;
    c)
	  config_file=$OPTARG
      ;;
    q)
	  quant_sw=$OPTARG
	    ;;
    t)
    technology=$OPTARG
      ;;
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done
# Paths
source $config_file



# array_length=($(wc -l ${docDir}samples_id.txt))

# Run quantification
cat ${docDir}samples_id.txt | while read in
do 
  #echo "qsub -P AG -A arubio -l h_vmem=6G -v config_file=$config_file,quant_sw=$quant_sw,sample_name=$in, technology=$technology -N '${in}_${quant_sw}_quantification' ${srcDir}02.quantification_united.sh"
  qsub -P AG -A arubio -l h_vmem=6G -v config_file=$config_file,quant_sw=$quant_sw,sample_name=$in,technology=$technology -N "${in}_${quant_sw}_quantification" ${srcDir}02.quantification_united.sh 
done
