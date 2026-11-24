#!/usr/bin/env bash
# Annotated analysis record; see README before running.
#Peak calling
# The original macs2.sh steps are shown below.
for file_name in *last.bam
do
bedtools bamtobed -i ${file_name} > ${file_name%bam}bed
macs2 callpeak -f BEDPE -g 119667750 -t ${file_name%bam}bed --nomodel --shift -100 --extsize 200 -n ${file_name%bam}last -q 0.05 --outdir ./peaks
done

# Calculate the percentage of input BED intervals intersecting peaks; verify interpretation.
# The original FRiP.sh steps are shown below.
ls *narrowPeak | while read id;
do
echo $id
bed=../$(basename $id "last_peaks.narrowPeak")bed
ls -lh $bed
Reads=$(bedtools intersect -a $bed -b $id | wc -l |awk '{print $1}')
totalReads=$(wc -l $bed|awk '{print $1}')
echo $Reads $totalReads
echo '==> FRiP value: ' $(bc <<< "scale=2;100*$Reads/$totalReads")'%'
done


