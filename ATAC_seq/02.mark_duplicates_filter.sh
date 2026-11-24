#!/usr/bin/env bash
# Annotated analysis record; see README before running.
# Mark duplicate alignments with samtools markdup; duplicates are not removed.
# The original markdup.sh steps are shown below.
#!/bin/bash

# Sort BAM files by read name.
for file_name in *_trimmed.bam
do 
   samtools sort -n -o ${file_name%_trimmed.bam}_dupsorted.bam ${file_name}
done

# Add mate information with fixmate.
for file_name in *_dupsorted.bam
do
   samtools fixmate -m ${file_name} ${file_name%dupsorted.bam}fixmated.bam
done

# Sort BAM files by genomic coordinate.
# Coordinate sorting follows.
for file_name in *_fixmated.bam
do
   samtools sort -o ${file_name%fixmated.bam}positionsort.bam ${file_name}
done

# Mark duplicate reads.
for file_name in *positionsort.bam
do
   samtools markdup ${file_name} ${file_name%positionsort.bam}markdup.bam
done
# End of original markdup.sh steps.

# Calculate alignment statistics; filter for proper pairs, MAPQ and organellar names.
ls *.bam | while read id;do(samtools index $id);done
ls *.bam | while read  id;do(samtools flagstat $id > ${id%bam}stat);done
ls *.bam | while read id;do(samtools view -@ 10 -h -f 2 -q 30 $id |grep -vE "Mt|Pt"|samtools view  -@ 10 -bS -o ${id%.rmdup.bam}.last.bam);done
ls *last.bam | while read id;do(samtools index $id);done


