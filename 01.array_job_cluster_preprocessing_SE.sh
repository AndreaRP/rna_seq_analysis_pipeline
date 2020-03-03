#!/bin/bash
set -e

#########################################################
#	         PREPROCESSING	
#########################################################
# 1. Performs quality trimming with Cutadapt.
# 2. Runs quality control of the trimmed reads with fastQC


# Config (paths)
source $config_file
# Read file
IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'

# Job id = algorithm.config (1 based)
jobid=${samples_id[$SGE_TASK_ID-1]}
sample_name=${jobid}


# FUNCTIONS
function makedir () {
	directory=$1
	if [ ! -d $directory ]
	then
		mkdir -p $directory
	fi
}
function filterReads (){
	# VARIABLES
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
	cutadapt -a AR00#=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR00#b=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AR000=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -a AR000b=GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -e 0.1 -O 5 -m 30 -o $samplePostProDir/${sample_name}_trimmed.fastq $rawDir/${sample_name}.fastq 2>&1 | tee -a $log
	echo -e "$(date) - Finished trimming sample ${sample_name}" >> $log
}

function filteredQC (){
	sample_name=$1
	# VARIABLES
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

# Call the functions
cmd1="filterReads"
cmd2="filteredQC"

#eval "$sample_name"
eval $cmd1
eval $cmd2
