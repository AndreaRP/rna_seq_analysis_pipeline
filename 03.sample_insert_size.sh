#!/bin/bash
set -e

#########################################################
#          INSERT SIZE METRICS     
#########################################################
# 1. Takes in job ID (sample name)
# 2. Runs bamtools for specified sample

# Config (paths)
source $config_file
# Read file (0 based)
IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'

# Job id = sample (1 based)
jobid=${samples_id[$SGE_TASK_ID-1]}
echo $jobid
# Access corresponding line in global file (0 based)

# Create the input (bam) and output file names 
bam_file=$(echo -e "${analysisDir}/04.quantification/RSEM/${jobid}/${jobid}.transcript.bam")
metrics_file=$(echo -e "${analysisDir}/05.insert_size/${jobid}_insert_size_metrics.txt")


# This part calls bamtools
cmd="bamtools stats -in ${bam_file} -insert > $metrics_file"

eval $cmd
