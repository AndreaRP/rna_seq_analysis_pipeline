#!/bin/bash
set -e

# Cluster version for benchmarking. Uses the algorithms_file.txt and configs_file.txt files in the src dir.
# Use: bash run_parallel.sh

while getopts "c:h" opt; do
  case $opt in
    h)
      echo "This will run the insert_size_array_job with the samples in doc/samples_id.txt.
            OPTIONS:
            -c: path to config file.
            Usage: bash ~/03.insert_size_array_job_submission.sh -c /path/to/config/file/file.config" 
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

if [ ! -d "${masterDir}sge_log" ]
  then
    mkdir -p "${masterDir}sge_log"
fi

if [ ! -d "${analysisDir}/05.insert_size" ]
  then
    mkdir -p "${analysisDir}/05.insert_size"
fi



# Use current working directory and current modules
#$ -cwd -V

IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'


# TASK_ID refers to each id in the arrayjob
qsub -P AG -A arubio -l h_vmem=5G -v config_file=$config_file -N "insert_size" -t 1:${#samples_id[@]} -tc 15 -e ${masterDir}sge_log/\$TASK_ID.e.log -o ${masterDir}/sge_log/\$TASK_ID.o.log ${srcDir}03.sample_insert_size.sh
