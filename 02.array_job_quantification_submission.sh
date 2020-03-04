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

IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'

qsub -P AG -A arubio -l h_vmem=5G -v config_file=$config_file,quant_sw=$quant_sw,technology=$technology -N "quantification" -t 1:${#samples_id[@]} -tc 10 -e ${analysisDir}/sge_log/\$TASK_ID.e.log -o ${analysisDir}/sge_log/\$TASK_ID.o.log ${srcDir}02.array_job_quantification.sh

