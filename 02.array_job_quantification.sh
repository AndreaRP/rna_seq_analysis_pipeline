#!/bin/bash
set -e

#########################################################
#	         QUANTIFICATION	
#########################################################
# Runs kallisto or rsem on specified samples by submitting 
# an array job.

# PARAMETERS	
# config_file: Configuration file. Typically docDir/parameters.config.
# Contains:
# 	masterDir: home directory for analysis
# 	analysisDir: dir containing processed files
# 	docDir: dir containing info related to analysis (config, sample list)
# 	rawDir: dir with raw data files (fastq)
# 	programsDir: dir with 3d party software
# 	srcDir: dir containing the scripts called in the analysis
# sample_name: Sample name (without extension)
# quant_sw: Quantification software to use (rsem|kallisto)
# technology: Sequencing technology (pe|se)


# Config (paths)
source $config_file
# Read file
IFS=$'\r\n' GLOBIGNORE='*' command eval  'samples_id=($(cat ${docDir}samples_id.txt))'

# Job id = algorithm.config (1 based)
jobid=${samples_id[$SGE_TASK_ID-1]}
sample_name=${jobid}
# Overwrite rawDir to trimmed reads directory
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
    # index_file="Mus_musculus.GRCm38.ercc96.84.chr.NoPseudoNomiRNA.rsem" # Normal option
    index_file="Mus_musculus.GCF_000001635.26_GRCm38.p6_genomic.rsem" # Just for liver bmal
    #index_file="Homo_sapiens.GRCh38.82.primary_assembly.NomiRNAPseudo.sorted.rsem"
    quantDir="${analysisDir}/04.quantification/RSEM/"
    # rsem_index="${refDir}/rsem/${index_file}" # Normal option
    rsem_index="/data3/genomes/mus_musculus/rsem/${index_file}" # Just for liver bmal
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
            echo -e "rsem-calculate-expression -p 8 --estimate-rspd --sort-bam-by-coordinate --output-genome-bam $fastq_file $rsem_index ${quantDir}/${sample_name}/${sample_name}" >> $log
            rsem-calculate-expression -p 8 --estimate-rspd --sort-bam-by-coordinate --output-genome-bam $fastq_file $rsem_index ${quantDir}/${sample_name}/${sample_name} 2>&1 | tee -a $log 
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

# Compose the commands
cmd1="makedir '${quantDir}/$sample_name/'"
cmd2="quantify"

# Call functions
eval $cmd1
eval $cmd2


