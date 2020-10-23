# RNA seq analysis (HPC)

Sequence of scripts to run necessary steps to obtain the counts in a bulk sequencing, for both SE and PE samples:
 - QC
 - Preprocesing
 - QC
 - Alinment
 - Quantification

The scripts run each step either in local machine, in an HPC as individual samples, or in an HPC as an arrayjob. 

## To run the pipeline

_In your own PC:_ 
### 1. Create the sample.txt list in DOC and the links to the raw fastq.gz files.
```bash
bash src/00.samples_gen_SE.sh -r ~project_dir//RAW/ -c ~project_dir//DOC/parameters.config
```
If you want to change the sample names, simply change the print parameters in line 71. Samples names cannot start with a number.

_In the cluster: (ssh user@nodename)_

### 2. Run the QC, trim with cutadapt and do new QC
```bash
bash /data3/user/src/01.cluster_preprocessing.sh -c ~project_dir//DOC/parameters.config -s se
```

### 3. Quantify with rsem / kallisto
```bash
bash /data3/user/src/02.cluster_quantification.sh -c ~project_dir//DOC/parameters.config -q rsem -t se
```

### 4. For PE, we can check the metrics by running
```bash
bash 03.insert_size_array_job_submission.sh -c ~project_dir/doc/parameters.config              
```
The machine number can be found in in the fastq headers:
@HWI-Mxxxx or @Mxxxx - MiSeq      
@HWUSI - GAIIx      
@HWI-Dxxxx - HiSeq 2000/2500      
@Kxxxx - HiSeq 3000(?)/4000      
@Nxxxx - NextSeq 500/550       
@Axxxxx - NovaSeq      