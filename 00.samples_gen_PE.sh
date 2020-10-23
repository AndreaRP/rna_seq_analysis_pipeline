#!/bin/bash
set -e 

#########################################################
#		  			SAMPLES GENERATOR				 	#
#########################################################

# Given raw read files with the following name structure: 
# sampleName_*R1*.fastq.gz, the sample name will be the word before the first '_' in the R1 reads file.
# 1. Gets the location of the parameters file
# 2. Gets the sample raw files location and creates symbolic link in 00-reads folder. 
# 3. Gets the sample name from the file names in 00-reads folder and creates samples_id.txt. 

# Input Files (In analysisDir/00-reads/)
# sampleName1_*R1*.fastq.gz
# sampleName2_*R1*.fastq.gz
# sampleName3_*R1*.fastq.gz
# ...

# Output Files (in analysisDir)
# samples_id.txt: File containing a list with the name of each sample.  

while getopts "c:r:h" opt; do
  case $opt in
    h)
      echo "This will generate the symbolic links to raw reads and the samples_id.txt.
OPTIONS:
  -c: path to config file.
  -r: path to raw files directory.
  -h: displays the help.
  Usage: bash ~/samples_gen.sh -r /path/to/raw/fastq/files/dir/ -c /path/to/config/file/" 
      ;;
    c)
    config_file=$OPTARG
      echo "Config file location: $config_file" 
      ;;
    r)
    raw_dir=$OPTARG
      echo "Raw reads directory location: $raw_dir" 
      ;;    
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

# CONSTANTS
source $config_file
readFolder="${analysisDir}00.reads/"

if [ ! -d ${readFolder} ]
then
  mkdir -p $readFolder
  echo -e "${readFolder} created"
fi

logFile="$readFolder/generation.log" 
echo -e "$(date) Generate Symbolic Links and samples_id.txt:" > $logFile
echo -e "Command: bash ~/samples_gen.sh -r $config_file -c $raw_dir" >> $logFile
echo -e "Config file location: $config_file" >> $logFile
echo -e "Raw reads directory location: $raw_dir" >> $logFile
echo -e "samples_id.txt: ${docDir}/samples_id.txt" >> $logFile

sampleName=""
> "${docDir}/samples_id.txt"
for file in $(find ${raw_dir}*RNA_Seq_*/ -name *.fastq.gz)
do
  sample=$(echo $file | rev | cut -d'/' -f1 | rev)

  # Generate samples_id.txt 
  sampleName=$(echo -e $sample | awk -F '_' -v OFS='_' '{print $1}')
  # Print only if it's not already been printed
  grep -q -F $sampleName "${docDir}/samples_id.txt" || echo $sampleName >> "${docDir}/samples_id.txt" 
  
  # Generate symbolic link
  sampleReadName=$(echo -e $sample | awk -F '_' -v OFS='_' '{print $1,$9}') # The link also has the R
  ln -s $file ${readFolder}/${sampleReadName}.fastq.gz
done

