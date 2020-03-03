
qualityDir="$analysisDir/03.trimmed_reads_QC/"
cd ${qualityDir}


# generate fastqc report: 
echo -e "Generate fastq report:" >> $lablog
echo -e "perl ./listFastQCReports.pl ${resultsDir}quality/data/ > ${resultsDir}quality/table.html" >> $lablog
#perl ./listFastQCReports.pl ./quality/data/ > ./quality/table.html
perl ./listFastQCReports.pl ${resultsDir}/data/ > ${qualityDir}/table.html
echo -e "perl ./createHTML.pl" >> $lablog
perl ./createHTML.pl