#!/usr/bin/env bash
# Annotated analysis record; see README before running.
# Generate RPKM-normalized bigWig tracks.
for file_name in *last.bam
do
samtools index -@ 10 ${file_name}
bamCoverage -b ${file_name} -of bigwig -o ../03.bw_file/${file_name%.bam}.bw -p 10 --binSize 10 --normalizeUsing RPKM
done


