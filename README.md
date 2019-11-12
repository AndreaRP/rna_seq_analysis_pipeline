# RNA seq analysis (HPC)

In your own PC: 
1. Create the sample.txt list in DOC and the links to the raw fastq.gz files.

bash src/00.samples_gen_SE.sh -r /data3/arubio/projects/Maria_WT_MerTK/RAW/ -c /data3/arubio/projects/Maria_WT_MerTK/DOC/parameters.config

If you want to change the sample names, simply change the print parameters in line 71. Samples names cannot start with a number.

In the cluster: (ssh arubio@samwise)

2. Run the QC, trim with cutadapt and do new QC

bash /data3/arubio/src/01.cluster_preprocessing.sh -c /data3/arubio/projects/Maria_WT_MerTK/DOC/parameters.config -s se

3. Quantify with rsem / kallisto

bash /data3/arubio/src/02.cluster_quantification.sh -c /data3/arubio/projects/Maria_WT_MerTK/DOC/parameters.config -q rsem -t se