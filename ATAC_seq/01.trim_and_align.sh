#!/usr/bin/env bash
# Annotated analysis record; see README before running.
# Trim adapters and low-quality reads; discard reads shorter than 20 bp.
ls *1.fastq.gz > 1
ls *2.fastq.gz > 2
paste 1 1 2 > config.raw

cat config.raw |while read id;
do echo $id
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
nohup trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired $fq1 $fq2 -o ../2.trimmed/ &
done

fastqc ./1.fastqc_trimmed/ *.fq.gz
ls *.fq.gz|while read id;do(gunzip $id);done

# Align reads to the reference genome with Bowtie2.
bowtie2-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.index
ls *1.fq > 1
ls *2.fq > 2
paste 1 1 2 > config.trimmed
cat config.trimmed | while read id; do echo $id
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
bowtie2 -x ./2.bowtie2_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.index -1 $fq1 -2 $fq2 -S ../4.bowtie2_file/${sample}.sam
done

# Convert SAM to BAM, sort and index.
ls *.sam | while read id;do(samtools view -bS $id > ${id%.sam}_trimmed.bam);done
ls *.bam|while read id;do(samtools sort -@ 8 $id -o ${id%trimmed.bam}sorted.bam);done
ls *_sorted.bam | while read id;do(samtools index -@ 10 $id);done


