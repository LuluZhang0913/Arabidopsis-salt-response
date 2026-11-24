# Annotated analysis record; see README before running.
# Differential accessibility analysis using DiffBind.
library(DiffBind)
samples <- read.table("./ATAC_sample_sheet.txt",head=T)
ATAC <- dba(sampleSheet=samples)
# Count reads over consensus peaks.
ATAC <- dba.count(ATAC, 
                 minOverlap=2,
                 bRemoveDuplicates=FALSE)
info <- dba.show(ATAC)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads*info$FRiP))
rownames(libsizes) <- info$ID
libsizes
ATAC <- dba.contrast(ATAC,categories=DBA_TREATMENT,minMembers = 2)
ATAC <- dba.analyze(ATAC,method=DBA_DESEQ2)
dba.show(ATAC,bContrasts=T)


