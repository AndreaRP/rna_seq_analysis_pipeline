#!/bin/bash
set -e

#########################################################
#	         PREPROCESSING SCRIPT			#
#########################################################
# 1. Creates necessary directories. 
# 2. Runs quality control of the raw reads with fastQC.
# 3. Performs quality trimming with Cutadapt.
# 4. Runs quality control of the trimmed reads with fastQC

# Parameters:
# -c config file
# -h = shows help

# Input Directory:
# rawDir= Directory where the raw reads are located (ANALYSIS/00.reads/)

# Input files: (In rawDir)
# sample_name.fastq.gz

# Output Directories:
# samplePreProQCDir = Raw reads quality analysis directory. (ANALYSIS/01.raw_reads_QC/sample_name)
# 	Output files: (In samplePreProQCDir)
#	sample_name_*_R1_*.fastqc.html: HTML file with the quality report for raw reads.
#	sample_name_*_R1_*.fastqc.zip: compressed file with the quality report files for raw reads generated by fastQC.
#	sample_name_preQC.log: log file for the quality analysis.
# samplePostProDir = Trimmed reads directory. (ANALYSIS/02.trimmed_reads/sample_name)
# 	Output files: (In samplePostProDir)
#	sample_name_R1_paired.fastq: fastq file with R1 trimmed paired reads.
#	sample_name_R1_unpaired.fastq: fastq file with R1 trimmed unpaired reads.
#	sample_name_postQC.log: log of the trimming process.
#	Note: Only the paired reads are used in the rest of the analysis
# samplePostProQCDir = Trimmed reads quality analysis directory. (ANALYSIS/03.trimmed_reads_QC/sample_name)
# 	Output files: (In samplePostProQCDir)
#	sample_name_R1_paired.fastqc.html: HTML file with the quality report for R1 trimmed reads.
#	sample_name_R1_paired.fastqc.zip: compressed file with the quality report files for R1 trimmed reads generated by fastQC.
#	sample_name_postQC.log: log file for the quality analysis.

# GET OPTS
while getopts "c:s:h" opt; do
  case $opt in
    h)
      echo "This will run FastQC on the raw reads.
OPTIONS:
	-c: path to config file.
	-h: displays the help.
	-s: name of the sample (as in samples_id.txt)
	Usage: bash ~/preprocessing.sh -c /path/to/config/file/file.config" 
      ;;
    c)
	  config_file=$OPTARG
      echo "Config file location: $config_file" 
      ;;	
    s)
	  sample_name=$OPTARG
      echo "Sample name: $sample_name" 
      ;;	
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

source $config_file
rawDir="$analysisDir/00.reads/"

#---------------------------------------------------------------------------------
# FUNCTIONS

function makedir () {
	directory=$1
	if [ ! -d $directory ]
	then
		mkdir -p $directory
	fi
}

function rawQC () {
	#	VARIABLES	
	samplePreProQCDir="$analysisDir/01.raw_reads_QC/$sample_name"
	log="$samplePreProQCDir/$sample_name_raw_QC.log"

	# CREATE DIRECTORY
	makedir $samplePreProQCDir

	# RAW READS QUALITY CONTROL
	> $log
	echo -e "$(date): Execute fastqc on $sample_name" >> $log
	echo -e " Command is: ### find $rawDir -name $sample_name.fastq* -exec fastqc {} --outdir $samplePreProQCDir \; ###" >> $log
	find $rawDir -name $sample_name.fastq* -exec fastqc {} --outdir $samplePreProQCDir \; >> $log
	echo -e "$(date): Finish fastqc" >> $log
	echo -e "$(date): Finished raw reads quaility control" >> $log
}
 
function filterReads (){
	#	VARIABLES
	samplePostProDir="$analysisDir/02.trimmed_reads/$sample_name"
	makedir $samplePostProDir	
	log="$samplePostProDir/${sample_name}_trimming.log"

	
	> $log	
	#-a: RNA-Seq TrueSeq y NEBNEXT adapter sequences
	#-e 0.1 : allow a mismatch rate of 1 mismatch in ten bases between the read and the adapters
	#-O 5 : The overlap must be at least 5 base-pairs
	#-m 30 : after trimming, reads less than 15bp are thrown out
	#-o sample_name_trimmed.fastq : put the trimmer data in this file
	echo -e "$(date) - Trimming sample ${sample_name}" >> $log
	echo -e "cutadapt -a AR00#=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR00#b=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR000=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -a AR000b=GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -e 0.1 -O 5 -m 30 -o ${sample_name}_trimmed.fastq $rawDir/${sample_name}.fastq.gz" >> $log
	cutadapt -a AR00#=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR00#b=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR000=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -a AR000b=GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -e 0.1 -O 5 -m 30 -o $samplePostProDir/${sample_name}_trimmed.fastq $rawDir/${sample_name}.fastq.gz 2>&1 | tee -a $log
	echo -e "$(date) - Finished trimming sample ${sample_name}" >> $log
}

function filteredQC (){
	sample_name=$1
	#	VARIABLES
	samplePostProQCDir="$analysisDir/03.trimmed_reads_QC/$sample_name"
	samplePostProDir="$analysisDir/02.trimmed_reads/$sample_name"
	makedir $samplePostProQCDir
	log="$samplePostProQCDir/${sample_name}_postQC.log"
	# TRIMMED READS QUALITY CONTROL
	echo -e "$(date): Execute fastqc.sh" >> "$log"
	find $samplePostProDir -name "*_trimmed.fastq" -exec fastqc {} --outdir $samplePostProQCDir \; >> "$log"
	echo -e "$(date): Finish fastqc.sh" >> "$log"
	echo -e "$(date): ********* Finished quaility control **********" >> "$log"
}

#---------------------------------------------------------------------------------
# RUN PREPROCESSING
filterReads 
filteredQC
# Create multiQC
multiqc $analysisDir/03.trimmed_reads_QC/*

