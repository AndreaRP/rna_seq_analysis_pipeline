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
while getopts "c:q:t:s:h" opt; do
  case $opt in
    h)
      echo "This will run Kallisto or rsem quantification on the processed reads.
OPTIONS:
  -h: displays the help.
	-c: path to config file.
	-s: Sample name.
  -q: Quantification software to use (rsem|kallisto)
  -t: Sequencing technology (pe|se)
	Usage: bash ~/quantification.sh -c /path/to/config/file/file.config -i index.file.name.idx" 
      ;;
    c)
	    config_file=$OPTARG
      echo "Config file location: $config_file" 
      ;;
    s)
    	sample_name=$OPTARG
    	echo "Sample Name: $sample_name"
	    ;;
    q)
      quant_sw=$OPTARG
      echo "Quantification: $quant_sw"
      ;;
    t)
      technology=$OPTARG
      echo "Technology: $technology"
      ;;
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

source $config_file
rawDir="${analysisDir}/02.trimmed_reads/"

function makedir () {
  directory=$1
  if [ ! -d "${directory}" ]
  then
    mkdir -p "${directory}"
  fi
}



case $quant_sw in
  rsem)
    index_file="Mus_musculus.GRCm38.ercc96.84.chr.NoPseudoNomiRNA.rsem"
    #index_file="Homo_sapiens.GRCh38.82.primary_assembly.NomiRNAPseudo.sorted.rsem"
    quantDir="${analysisDir}/04.quantification/RSEM/"
    rsem_index="${refDir}/rsem/${index_file}"
    ;;
  kallisto)
    # index_file="Mus_musculus.GRCm38.ercc96.84.chr.NoPseudoNomiRNA.kallisto.idx"
    index_file="mm_GRCm38_noPseudoNomiRNA_GFP.kallisto.idx"
    quantDir="${analysisDir}/04.quantification/kallisto/"
    kallisto_index="${refDir}/kallisto/${index_file}"
    ;;
  :)
    echo "Option $OPTARG requires an argument." 
    exit 1
    ;;
esac

log="${quantDir}/${sample_name}/${sample_name}.log"

# function create_index {
#   # Create index file
#   if [ ! -a "${refDir}/${index_file}.idx" ]
#     then
#       ${programsDir}/kallisto index -i "${refDir}/${index_file}.idx" ${refDir}/$index_file  2>&1 | tee -a "${docDir}/{$index_file}.log"
#     fi
# }

function quantify () {
case $quant_sw in
  kallisto)
      case $technology in 
        se)            
            input_file="${rawDir}/${sample_name}/${sample_name}_trimmed.fastq"            
            # Run kallisto with options:
            #-i: $index_file.idx
            #-o: "${analysisDir}/04-quantification/sample_name/" 
            #--bootstrap-samples=100 (For later Diff expression analysis)
            #--single: single-end reads
            #--fragment-length=200 
            #--sd=50 
            echo -e "$(date) - $sample_name" > $log
            echo -e "${programsDir}/kallisto quant -i '${kallisto_index}' -o '${quantDir}/${sample_name}/' -b 100 --single --fragment-length 200 --sd 50 $input_file" >> $log
            ${programsDir}/kallisto quant -i "${kallisto_index}" -o "${quantDir}/${sample_name}/" -b 100 --single --fragment-length 200 --sd 50 $input_file
          #2>&1 | tee -a $log &
            echo -e "DONE" >> $log
        ;;
        pe)
            input_file_1="${rawDir}/${sample_name}/${sample_name}_R1_trimmed.fastq"
            input_file_2="${rawDir}/${sample_name}/${sample_name}_R2_trimmed.fastq"          
            # Run kallisto with options:
            #-i: $index_file.idx
            #-o: "${analysisDir}/04-quantification/sample_name/" 
            #--bootstrap-samples=100 (For later Diff expression analysis)
            #--fragment-length=200 
            #--sd=50 
            echo -e "$(date) - $sample_name" > $log
            echo -e "${programsDir}/kallisto quant -i '${kallisto_index}' -o '${quantDir}/${sample_name}/' -b 100 --single --fragment-length 200 --sd 50 $input_file_1 $input_file_2" >> $log
            ${programsDir}/kallisto quant -i "${kallisto_index}" -o "${quantDir}/${sample_name}/" -b 100 --fragment-length 200 --sd 50 $input_file_1 $input_file_2
          #2>&1 | tee -a $log &
            echo -e "DONE" >> $log
        ;;
      esac
  ;;
  rsem)
      case $technology in 
        se)
            fastq_file="${rawDir}/${sample_name}/${sample_name}_trimmed.fastq"
            echo -e "$(date) - $sample_name" > $log
            echo -e "rsem-calculate-expression -p 8 --estimate-rspd --output-genome-bam $fastq_file $rsem_index ${quantDir}/${sample_name}/${sample_name}" >> $log
            rsem-calculate-expression -p 8 --estimate-rspd --output-genome-bam $fastq_file $rsem_index ${quantDir}/${sample_name}/${sample_name} 2>&1 | tee -a $log 
            echo -e "DONE" >> $log
        ;;
        pe)
            input_file_1="${rawDir}/${sample_name}/${sample_name}_R1_trimmed.fastq"
            input_file_2="${rawDir}/${sample_name}/${sample_name}_R2_trimmed.fastq" 
            echo -e "$(date) - $sample_name" > $log
            echo -e "rsem-calculate-expression -p 8 --estimate-rspd --output-genome-bam $fastq_file $rsem_index ${quantDir}/${sample_name}/${sample_name}" >> $log
            rsem-calculate-expression -p 8 --estimate-rspd --output-genome-bam --paired-end $input_file_1 $input_file_2 $rsem_index ${quantDir}/${sample_name}/${sample_name} 2>&1 | tee -a $log 
            echo -e "DONE" >> $log
        ;;
        :)
            echo "Option $OPTARG requires an argument." 
            exit 1
        ;;    
      esac
  ;;
  :)
    echo "Option $OPTARG requires an argument." 
    exit 1
  ;;
esac
}

# Run the script
	makedir "${quantDir}/$sample_name/"    
	quantify


