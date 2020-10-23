#!/bin/bash
set -e 

#########################################################
#           SAMPLES GENERATOR         #
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

while getopts "r:s:f:h" opt; do
  case $opt in
    h)
      echo "This will generate the symbolic links to raw reads and the metadata file
OPTIONS:
  -s: path to submission dir.
  -r: path to raw files directory.
  -f: name of the metadata file.
  -h: displays the help.
  Usage: bash ~/samples_gen.sh -r /path/to/raw/fastq/files/dir/ -s /path/to/submission/dir/ -f metadata" 
      ;;
    r)
    raw_dir=$OPTARG
      echo "Raw reads directory location: $raw_dir" 
      ;;  
    s)
    sub_dir=$OPTARG
      echo "Directory to write the link to (SUBMISSION): $sub_dir" 
      ;;  
    f)
    metadata_info=$OPTARG
      echo "Name for the metadata file: $metadata_info" 
      ;;  
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done


if [ ! -d ${sub_dir} ]
then
  mkdir -p $sub_dir
  echo -e "${sub_dir} created"
fi


sampleName=""
> "${sub_dir}/${metadata_info}.csv"
for file in $(find ${raw_dir}* -name *.fastq.gz)
do
  # Get sample name
  sampleName=$(echo $file | rev | cut -d'/' -f1 | rev)
  sampleName=$(echo -e $sampleName | awk -F '_' -v OFS='_' '{print $1}')

  # Generate symbolic link to submission dir
  fileName=$(basename ${file})
  ln -s $file "${sub_dir}/${fileName}"

    # METADATA
  # Get checksum
  md5=$(md5sum $file | cut -f1 -d' ')
  echo "$sampleName,$fileName,$md5" >> "${sub_dir}/${metadata_info}.csv"
done

