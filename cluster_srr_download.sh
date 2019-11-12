#!/bin/bash
set -e

while getopts "i:o:h" opt; do
  case $opt in
    h)
      echo "This will download the requested files.
OPTIONS:
	-i: path to srr list
	-o: path to output location
	-h: Runs help
	Usage: bash ~/cluster_srr_download.sh -i /path/to/srr/ftp/url/list.txt -o /path/to/output/dir/" 
	exit 1
      ;;
    i)
	url_list=$OPTARG
      ;;
    o)
	output=$OPTARG
	  ;;
    :)
      echo "Option $OPTARG requires an argument." 
      exit 1
      ;;
  esac
done

n_jobs=$(wc -l $url_list)

  cat $url_list | awk '{ print $10 }' | while read in
do
	qsub -P AG -A arubio -l h_vmem=2G -N "fastq-dump_${in}" /data3/arubio/src/srr_download.sh -z ${in} -o ${output} -t 0-${n_jobs}%2
done


