# example using bedtools genomecov + RADAR to find out number of rna editings sites with coverage
options(java.parameters = "-Xmx20g")

library(rtracklayer)
library(stringr)
library(plyr)
library(dplyr)
library(xlsx)

radarFile <- "Human_AG_all_hg19_v2.txt"
darnedFile <- "DARNED.txt"

radar <- read.table(radarFile, header=T)
darned <- read.table(darnedFile, header=F, skip=1, sep="\t")
names(darned) <- c("Chrom",	"Coordinate", "Strand", "DNA", "RNA",	"Gene",	"Location", "ExReg", "Source", "PubMed ID")

# fix contig names
radar$chromosome <- sub("chr","",radar$chromosome)
radar$chromosome[radar$chromosome=="M"] <- "MT" 
radarRanges <- with(radar, GRanges(seq=chromosome, ranges=IRanges(start=position,end=position), strand=strand, gene=gene,
                                   annot1=annot1, annot2=annot2, alu=alu., non_alu_repetitive=non_alu_repetitive.,
                                   conservation_chimp=conservation_chimp, conservation_rhesus,conservation_rhesus,
                                   conservation_mouse=conservation_mouse))

# fix the strand stuff with darned
darned$Strand <- sub("?","*", darned$Strand, fixed=T)
darnedRanges <- with(darned, GRanges(seq=Chrom, ranges=IRanges(start=Coordinate, end=Coordinate), strand=Strand, dna=DNA, rna=RNA,
                                     gene=Gene, location=Location, exreg=ExReg, source=Source, pubmed="PubMed ID"))

# sort stuff so the seqlevels actually match what they should
radarRanges <- sortSeqlevels(radarRanges)
radarRanges <- sort(radarRanges)
darnedRanges <- sortSeqlevels(darnedRanges)
darnedRanges <- sort(darnedRanges)

# remove irrational coordinates
darnedRanges <- darnedRanges[!(seqnames(darnedRanges)=="17" & start(darnedRanges) > 83257441)]

# example case
# coverage calculated w/ bedtools genomecov -bga -split
#bedgraph <- import.bedGraph("rna_01pg_F-P21-M1-65-1_run72_016_Aligned.sortedByCoord.out.bedgraph.gz")

# mcols(bedgraph[findOverlaps(radarRanges, bedgraph, select="first"),])$score
#shortname <- str_match(foo,"rna_(.*)_Aligned.*")[2]
#radarRanges$shortname <- mcols(bedgraph[findOverlaps(radarRanges, bedgraph, select='first'),])$score
#mcols(radarRanges) <- rename(mcols(radarRanges), c(shortname=shortname))
#darnedRanges$shortname <- mcols(bedgraph[findOverlaps(darnedRanges,bedgraph, select="first"),])$score
#mcols(darnedRanges) <- rename(mcols(darnedRanges), c(shortname=shortname))

# we only want the aligned files
files <- dir(path="/data/projects/rna_prep_comparison/subsampled_10M/star", pattern=".*_Aligned.*.bedgraph.gz",full.names = T)
shortnames <- sapply(files, function(x) str_match(basename(x),"rna_(.*)_Aligned.*")[2])

# yes, this is stupid, to transpose, but there is something bizarre memory wise when I use lapply derivatives
tbl <- ldply(files, function(x) {
  bedgraph <- import.bedGraph(x)
  mcols(bedgraph[findOverlaps(radarRanges, bedgraph, select='first'),])$score
})
rownames(tbl) <- shortnames
tbl <- as.data.frame(t(tbl))
mcols(radarRanges) <- cbind(mcols(radarRanges), tbl)

dtbl <- ldply(files, function(x) {
  bedgraph <- import.bedGraph(x)
  mcols(bedgraph[findOverlaps(darnedRanges, bedgraph, select='first'),])$score
})
rownames(dtbl) <- shortnames
dtbl <- as.data.frame(t(dtbl))
mcols(darnedRanges) <- cbind(mcols(darnedRanges), dtbl)

# remove all zero rows to get the table size down
# too many rows pass. no point

# 5 samples 5 reads
filteredRadarRanges <- radarRanges[apply(tbl,1,function(x) sum(x>=5)>=5),]
filteredDarnedRanges <- darnedRanges[apply(dtbl,1,function(x) sum(x>=5)>=5),]
easytbl <- tbl[apply(tbl,1,function(x) sum(x>=5)>=5),]
easydtbl <- dtbl[apply(dtbl,1,function(x) sum(x>=5)>=5),]

# count up sites with more than 5 reads
summarydf <- data.frame(
  "radar"=apply(easytbl, 2, function(x) sum(x>=5)),
  "darned"=apply(easydtbl, 2, function(x) sum(x>=5)))

write.csv(summarydf, "editingSiteCoverage_10M.csv")
write.csv(as.data.frame(filteredRadarRanges), "radarSiteCoverage_10M.csv")
write.csv(as.data.frame(filteredDarnedRanges), "darnedSiteCoverage_10M.csv")

# xlsx is taking a geologic amount of time, just combine the csv's
write.xlsx(summarydf,"editingSiteCoverage.xlsx", sheetName = "Summary")
write.xlsx(as.data.frame(filteredRadarRanges), "editingSiteCoverage.xlsx", sheetName = "Radar sites", append=T)
write.xlsx(as.data.frame(filteredDarnedRanges), "editingSiteCoverage.xlsx", sheetName = "Darned sites", append=T)