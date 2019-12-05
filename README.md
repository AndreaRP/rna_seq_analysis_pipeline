# RNA seq analysis (HPC)

In your own PC: 
1. Create the sample.txt list in DOC and the links to the raw fastq.gz files.

bash src/00.samples_gen_SE.sh -r ~project_dir//RAW/ -c ~project_dir//DOC/parameters.config

If you want to change the sample names, simply change the print parameters in line 71. Samples names cannot start with a number.

In the cluster: (ssh user@samwise)

2. Run the QC, trim with cutadapt and do new QC

bash /data3/user/src/01.cluster_preprocessing.sh -c ~project_dir//DOC/parameters.config -s se

3. Quantify with rsem / kallisto

bash /data3/user/src/02.cluster_quantification.sh -c ~project_dir//DOC/parameters.config -q rsem -t se


The machine number can be found in in the fastq headers:

@HWI-Mxxxx or @Mxxxx - MiSeq      
@HWUSI - GAIIx      
@HWI-Dxxxx - HiSeq 2000/2500      
@Kxxxx - HiSeq 3000(?)/4000      
@Nxxxx - NextSeq 500/550       
@Axxxxx - NovaSeq      