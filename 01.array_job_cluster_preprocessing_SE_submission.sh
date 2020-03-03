#!/bin/bash
set -e

# Cluster version for benchmarking. Uses the algorithms_file.txt and configs_file.txt files in the src dir.
# Use: bash run_parallel.sh

while getopts "c:h" opt; do
  case $opt in
    h)
      echo "This will run 01.array_job_cluster_preprocessing with the samples in doc/samples_id.txt.
            OPTIONS:
            -c: path to config file.
            Usage: bash ~/01.array_job_cluster_preprocessing_submission.sh -c /path/to/config/file/file.config" 
      ;;
    c)
	config_file=$OPTARG
      ;;
    :)
    	echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

source $config_file

if [ ! -d "${analysisDir}/sge_log" ]
  then
    mkdir -p "${analysisDir}/sge_log"
fi



# Use current working directory and current modules
#$ -cwd -V

IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'


# TASK_ID refers to each id in the arrayjob (tc=maximum number of simultaneous jobs)
qsub -P AG -A arubio -l h_vmem=5G -v config_file=$config_file -N "preprocessing" -t 1:${#samples_id[@]} -tc 5 -e ${analysisDir}/sge_log/\$TASK_ID.e.log -o ${analysisDir}/sge_log/\$TASK_ID.o.log ${srcDir}01.array_job_cluster_preprocessing_SE.sh
#echo "qsub -P AG -A arubio -l h_vmem=5G -v config_file=$config_file -N "preprocessing" -t 1:${#samples_id[@]} -tc 5 -e ${analysisDir}/sge_log/\$TASK_ID.e.log -o ${analysisDir}/sge_log/\$TASK_ID.o.log ${srcDir}01.array_job_cluster_preprocessing_SE.sh"
# Create multiQC once all QC is done
# multiqc $analysisDir/03.trimmed_reads_QC/*
