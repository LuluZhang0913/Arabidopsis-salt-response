# Annotated analysis record; see README before running.
##ChIPseeker 
library(ChIPseeker)
ath <- txdbmaker::makeTxDbFromGFF('Arabidopsis_thaliana.TAIR10.54.gff3')
peak <- readPeakFile("all_peaks.bed", skip = 1)
length(peak)
peakAnno <- annotatePeak(peak, tssRegion = c(-2000, -1), TxDb = ath)
write.table(peakAnno, file = "anno_all_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE)
plotAnnoBar(peakAnno)
plotAnnoPie(peakAnno)

