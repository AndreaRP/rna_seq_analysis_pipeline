#!/bin/bash
set -e

#########################################################
#		  			PREPROCESSING SCRIPT			 	#
#########################################################
# 1. Creates necessary directories. 
# 2. Creates index files if necessary
# 3. Runs kallisto to quantify.

# Parameters:
# -c config file
# -i index file
# -s sample name
# -h = shows help

# Input Directory:
# rawDir= Directory where the raw reads are located (ANALYSIS/00-reads/)

# Input files: (In rawDir)
# sample_name.fastq.gz

# Output Directories:
# quantDir = Raw reads quality analysis directory. (ANALYSIS/04-quantification/sample_name)
# 	Output files: (In QuantDir)
#	sample_name_abundance.txt: Abundances are reported in “estimated counts” (est_counts) and in Transcripts Per Million (TPM) in tabbed format.


# GET OPTS
while getopts "c:i:s:h" opt; do
  case $opt in
    h)
      echo "This will run Kallisto quantification on the processed reads.
OPTIONS:
	-c: path to config file.
	-h: displays the help.
  	-i: Name of the transcriptome file. Must be in $refDir
	-s: Sample name.
	Usage: bash ~/quantification.sh -c /path/to/config/file/file.config -i index.file.name.idx" 
      ;;
    c)
	config_file=$OPTARG
      	echo "Config file location: $config_file" 
      ;;
    i)
    	index_file=$OPTARG
      	echo "Index file: $index_file" 
      ;;
    s)
	sample_name=$OPTARG
	echo "Sample Name: $sample_name"
	;;
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

source $config_file
rawDir="${analysisDir}/02.trimmed_reads/"
quantDir="${analysisDir}/04.quantification/kallisto/"


function makedir () {
  directory=$1
  if [ ! -d "${directory}" ]
  then
    mkdir -p "${directory}"
  fi
}

# function create_index {
#   # Create index file
#   if [ ! -a "${refDir}/${index_file}.idx" ]
#     then
#       ${programsDir}/kallisto index -i "${refDir}/${index_file}.idx" ${refDir}/$index_file  2>&1 | tee -a "${docDir}/{$index_file}.log"
#     fi
# }

function quantify () {
  #sample_name=$1
  log="${quantDir}/${sample_name}/${sample_name}.log"
  output_file="${rawDir}/${sample_name}/${sample_name}_trimmed.fastq"
  kallisto_index="${refDir}/kallisto/${index_file}"
  # Run kallisto with options:
  #-i: $index_file.idx
  #-o: "${analysisDir}/04-quantification/sample_name/" 
  #--bootstrap-samples=100 (For later Diff expression analysis)
  #--single: single-end reads
  #--fragment-length=200 
  #--sd=50 
  echo -e "$(date) - $sample_name" > $log
  echo -e "${programsDir}/kallisto quant -i '${kallisto_index}' -o '${quantDir}/${sample_name}/' -b 100 --single --fragment-length 200 --sd 50 $output_file" >> $log
  ${programsDir}/kallisto quant -i "${kallisto_index}" -o "${quantDir}/${sample_name}/" -b 100 --single --fragment-length 200 --sd 50 $output_file
#2>&1 | tee -a $log &
  echo -e "DONE" >> $log
}

# Run the script
	makedir "${quantDir}/$sample_name/"    
	quantify


