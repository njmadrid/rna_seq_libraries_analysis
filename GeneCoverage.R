### GeneCoverage.R
###
###    Genome-wide exon coverage per gene using RNA-seq data.
###
### Copyright 2017 Nathaniel Madrid
###
###
### This file is made available under the terms of the GNU General Public License, 
### version 3, or at your option, any later version, incorporated herein by reference.

library(DBI); # Required for Homo.sapiens
library(Homo.sapiens);
library(plyr);
library(Rsamtools);

NUMBINS = 100;
NUMCORES = 12;

HtseqDir = "../../RNAPrepSamples/htseq_exons_tmap"; # htseq files used to quickly identify highly expressed genes.
HtseqLocation = list.files(path=HtseqDir, full.names=T);
HtseqLocation = "../../../../projects/rna_prep_comparison/htseq_exons_tmap/rna_10ng_SMARTer_run49_004_feature_counts.txt";
BamLocation = "../../../../projects/rna_prep_comparison/star/rna_10ng_SMARTer_run49_004_Aligned.sortedByCoord.out.bam";

htseq = read.table(file=HtseqLocation, header=T);
htseq$Gene = gsub("_exon_[ -~]+", "", htseq$feature_id);
FeatureTable = aggregate(htseq$feature_id, by=list(Gene=htseq$Gene), FUN=length);
FeatureTable = cbind(FeatureTable[ , 2], aggregate(htseq$RPKM, by=list(Gene=htseq$Gene), FUN=max));
colnames(FeatureTable) <- c("NumExons", "Gene", "MaxRPKM");
FeatureTable$symbol = (select(Homo.sapiens, FeatureTable$Gene, "SYMBOL", "REFSEQ"))$SYMBOL;
FeatureTable$GeneID = (select(Homo.sapiens, FeatureTable$Gene, "GENEID", "REFSEQ"))$GENEID;
Strandedness = select(Homo.sapiens, htseq$Gene, "EXONSTRAND", "REFSEQ");
# Ignore warning. For each gene, a strand value is returned for per exon ("GENESTRAND" 
# is not an available column). The ddply command will produce one strand value per gene.
Strandedness$EXONSTRAND[is.na(Strandedness$EXONSTRAND)] <- "*";

GetGeneStrand = function(ExonStrands) {
  ExonStrandTable = table(ExonStrands);
  ExonStrandTableName = names(table(ExonStrands));
  ExonStrandTableValue = table(ExonStrands);
  return(ExonStrandTableName[which.max(ExonStrandTableValue)]);
}

FeatureTable$Strand = (ddply(Strandedness, .(REFSEQ), summarize, strand=GetGeneStrand(EXONSTRAND)))$strand;
FinalTable = FeatureTable[FeatureTable$MaxRPKM >= 50 & FeatureTable$Strand %in% c("+", "-"), ]; # MaxRPKM = 50 is too small

AllExons = read.table(file="RefSeq_hg19_exons_nodups_021114g1k.bed", header=T);
colnames(AllExons) <- c("chr", "start", "end", "refSeq", "hg19", "strand", "description", "other");
FinalTable$Gene = paste(FinalTable$Gene, "_", sep="");
TopExons = mclapply(FinalTable$Gene, function(g) grep(g, AllExons$refSeq), mc.cores=NUMCORES);
TopExons = mclapply(TopExons, function(g) data.frame(chr=AllExons$chr[g], start=AllExons$start[g], end=AllExons$end[g]), mc.cores=NUMCORES);
TopGR = mclapply(TopExons, function(g) GRanges(seqnames=g$chr, ranges=IRanges(start=g$start, end=g$end)), mc.cores=NUMCORES);

getPileUpPerGene = function(which) {
  pu = pileup(BamLocation,
              scanBamParam=ScanBamParam(which=which),
              pileupParam=PileupParam(max_depth=10000,
                                      distinguish_strands=F,
                                      distinguish_nucleotides=F,
                                      include_deletions=F,
                                      include_insertions=F));
  return(pu);
}

PileUpPerGene = mclapply(TopGR, getPileUpPerGene, mc.cores=NUMCORES);
HasPileUps = unlist(mclapply(PileUpPerGene, function(d) nrow(d)>0, mc.cores=NUMCORES));
puUseful = PileUpPerGene[HasPileUps];
UsefulTable = FinalTable[HasPileUps, ];

MergeByExon = function(PileIndex) {
  p = puUseful[[PileIndex]];
  MergedDF = data.frame(pos=p$pos, count=p$count);
  MergedDF[is.na(MergedDF)] <- 0;
  MergedDF = MergedDF[order(MergedDF$pos), ];
  return(MergedDF);
}

puUnfolded = mclapply(1:length(puUseful), MergeByExon, mc.cores=NUMCORES);

isValidCount = function(pos, count) {
  tryCatch(rep(pos, count),
           warning = function(w) { return(F) },
           error = function(e) { return(F) }) 
  return(T);
}

BinByExon = function(i) {
  d = puUnfolded[[i]];
  ValidCounts = lapply(1:nrow(d), function(j) isValidCount(d$pos[j], d$count[j]));
  d = d[isValidCount(d)];
  d$pos = c(1:length(d$pos)) / length(d$pos) * NUMBINS;
  breaks = 0:NUMBINS;
  flattened = unlist(lapply(1:nrow(d), function(j) rep(d$pos[j], d$count[j])));
  HistCounts = hist(flattened, breaks=breaks, right=F, plot=F)$counts;
  return(HistCounts);
}

puBin = mclapply(1:length(puUnfolded), BinByExon, mc.cores=NUMCORES);
directionalBin = mclapply(1:length(puBin), function(i) { if(UsefulTable$Strand[i]=="-") { return(rev(unlist(puBin[[i]]))) } else { return(puBin[[i]]) }}, mc.cores=NUMCORES);
normBin = mclapply(directionalBin, function(x) scale(x), mc.cores=NUMCORES);

BinFlat = data.frame(count=unlist(normBin));
BinFlat$BinID = rep(1:NUMBINS, length.out=nrow(BinFlat));
BinFlat$Gene = rep(UsefulTable$Gene, rep(NUMBINS * length(BamLocation), nrow(UsefulTable)));
BinFlat$Method = gsub(paste(HtseqDir, "/|rna_|_feature_counts.txt", "../../../../projects/rna_prep_comparison/htseq_exons_tmap/rna_", sep="|"), "", HtseqLocation);

# For oligo dT methods only!
ExpandedBinValues = mclapply(directionalBin, function(x) unlist(lapply(1:NUMBINS, function(i) rep(i, unlist(x)[i]))), mc.cores=NUMCORES);
centroids = unlist(mclapply(ExpandedBinValues, mean, mc.cores=NUMCORES));
AcceptableGenes = UsefulTable$Gene[which(centroids>50)];
BinFlat = BinFlat[which(BinFlat$Gene %in% AcceptableGenes), ];

# For each bin of each sample, get the mean count across all genes.
BinPerFile = ddply(BinFlat, .(BinID, Method), summarise, zscore=mean(count)); # mean
ScaledMinZero = BinPerFile$zscore + abs(min(BinPerFile$zscore));
BinPerFile$RescaledCount = ScaledMinZero / max(ScaledMinZero);
BinPerFile$Method = factor(BinPerFile$Method, labels=BinFlat$Method[1]);

g = ggplot(BinPerFile, aes(x=BinID, y=RescaledCount, color=Method)) + geom_line();
g = g + ylab(label="Normalized Read Coverage") + xlab("Normalized Distance Along Transcript 5' -> 3' (%)");
g = g + ggtitle(paste(format(length(unique(BinFlat$Gene)), big.mark=","), " Genes", sep=""));
#g;

####################################################################################
# Continue after an appropriate set of genes (i.e. some minimum number of reads) has 
# been identified. These genes will then be used to create the gene coverage plot(s).
####################################################################################

library(GenomicRanges);
library(Rsamtools);
library(plyr);
library(pheatmap);
library(ggplot2);
library(data.table);
library(parallel);

BamDirectory = "../../../RNAPrepSamplesFINAL/nodups";
FileLocations = list.files(path=BamDirectory, pattern="\\.bam$", full.names=T);

#SampleNames = c("SMARTer Oligo dT", "MALBAC", "CD34 Custom 1", "CD34 Custom 2", "CD34 Life Tech 1", "CD34 Life Tech 2", "Life Tech 1", "Custom 1", "Custom 2", "Life Tech 2");
#SampleOrdering = c(3:6, 8, 9, 7, 10, 2, 1);
FinalTable = data.frame(Gene=unique(BinFlat$Gene));
FinalTable$symbol = (select(Homo.sapiens, sub("_$", "", as.character(FinalTable$Gene)), "SYMBOL", "REFSEQ"))$SYMBOL;
FinalTable$GeneID = (select(Homo.sapiens, sub("_$", "", as.character(FinalTable$Gene)), "GENEID", "REFSEQ"))$GENEID;
Strandedness = select(Homo.sapiens, sub("_$", "", as.character(FinalTable$Gene)), "EXONSTRAND", "REFSEQ");
# Ignore warning. For each gene, a strand value is returned for per exon ("GENESTRAND" 
# is not an available column). The ddply command will produce one strand value per gene.
Strandedness$EXONSTRAND[is.na(Strandedness$EXONSTRAND)] <- "*";

GetGeneStrand = function(ExonStrands) {
  ExonStrandTable = table(ExonStrands);
  ExonStrandTableName = names(table(ExonStrands));
  ExonStrandTableValue = table(ExonStrands);
  return(ExonStrandTableName[which.max(ExonStrandTableValue)]);
}

FinalTable$Strand = (ddply(Strandedness, .(REFSEQ), summarize, strand=GetGeneStrand(EXONSTRAND)))$strand;

AllExons = read.table(file="../RefSeq_hg19_exons_nodups_021114g1k.bed", header=T);
colnames(AllExons) <- c("chr", "start", "end", "refSeq", "hg19", "strand", "description", "other");
TopExons = mclapply(FinalTable$Gene, function(g) grep(g, AllExons$refSeq), mc.cores=12);
TopExons = mclapply(TopExons, function(g) data.frame(chr=AllExons$chr[g], start=AllExons$start[g], end=AllExons$end[g]), mc.cores=12);
TopGR = mclapply(TopExons, function(g) GRanges(seqnames=g$chr, ranges=IRanges(start=g$start, end=g$end)), mc.cores=12);

getPileUpPerGene = function(which) {
  pu = lapply(FileLocations,
              function(x) pileup(x,
                                 scanBamParam=ScanBamParam(which=which),
                                 pileupParam=PileupParam(max_depth=10000,
                                                         distinguish_strands=F,
                                                         distinguish_nucleotides=F,
                                                         include_deletions=F,
                                                         include_insertions=F)));
  return(pu);
}

# ~5-30 minutes
date()
PileUpPerGene = mclapply(TopGR, getPileUpPerGene, mc.cores=12);
date()

MergeByExon = function(PileIndex) {
  NumFiles = length(PileUpPerGene[[PileIndex]]);
  p = PileUpPerGene[[PileIndex]][[1]];
  MergedDF = data.frame(pos=p$pos, count=p$count);
  for(i in 2:NumFiles) {
    p = PileUpPerGene[[PileIndex]][[i]];
    MergedDF = merge(MergedDF, data.frame(pos=p$pos, count=p$count), by="pos", all=T);
  }
  colnames(MergedDF)[seq(from=2, by=1, length.out=NumFiles)] <- c(paste("count", c(1:NumFiles), sep=""));
  MergedDF[is.na(MergedDF)] <- 0;
  MergedDF = MergedDF[order(MergedDF$pos), ];
  return(MergedDF);
}

puUnfolded = mclapply(1:length(PileUpPerGene), MergeByExon, mc.cores=12);

isValidCount = function(pos, count) {
  tryCatch(rep(pos, count),
           warning = function(w) { return(F) },
           error = function(e) { return(F) }) 
  return(T);
}

BinByExon = function(i) {
  d = puUnfolded[[i]];
  NumFiles = length(d);
  HistList = list();
  for(f in 2:NumFiles) {
    ValidCounts = lapply(1:nrow(d), function(j) isValidCount(d$pos[j], d[ , f][j]));
    d = d[isValidCount(d)];
    d$pos = c(1:length(d$pos)) / length(d$pos) * NUMBINS;
    breaks = 0:NUMBINS;
    flattened = unlist(lapply(1:nrow(d), function(j) rep(d$pos[j], d[ , f][j])));
    HistCounts = hist(flattened, breaks=breaks, right=F, plot=F)$counts;
    HistList[[f-1]] = HistCounts;
  }
  return(HistList);
}

puBin = mclapply(1:length(puUnfolded), BinByExon, mc.cores=12);

ReverseOrder = function(i) { 
  FilesPerGene = puBin[[i]];
  NumFiles = length(FilesPerGene);
  ReturnList = list();
  if(FinalTable$Strand[i]=="-") { 
    for(f in 1:NumFiles) {
      ReturnList[[f]] = rev(unlist(FilesPerGene[f])); 
    }
  } 
  else { 
    for(f in 1:NumFiles) {
      ReturnList[[f]] = unlist(FilesPerGene[f]);
    }
  }
  return(ReturnList);
}

directionalBin = mclapply(1:length(puBin), ReverseOrder, mc.cores=12);
normBin = mclapply(directionalBin, function(x1) lapply(x1, function(x2) scale(x2)), mc.cores=12);

BinFlatV2 = data.frame(count=unlist(normBin));
BinFlatV2$BinID = rep(1:NUMBINS, length.out=nrow(BinFlatV2));
FileNames = gsub("../../../RNAPrepSamplesFINAL/nodups/rna_|_Aligned.sortedByCoord.out.bam", "", FileLocations);
BinFlatV2$Gene = rep(FinalTable$Gene, rep(NUMBINS * length(FileNames), nrow(FinalTable)));
BinFlatV2$Method = rep(rep(FileNames, rep(NUMBINS, length(FileNames))), length(normBin));

# For each bin of each sample, get the mean count across all genes.
MeanWithNA = function(counts) { counts[is.na(counts)]<-0; return(mean(counts)) } # mean with NAs
BinPerFile = ddply(BinFlatV2, .(BinID, Method), summarise, zscore=MeanWithNA(count)); # mean with NAs

ScaleMinMax = function(bin) {
  ScaledMinZero = bin$zscore + abs(min(bin$zscore));
  bin$RescaledCount = ScaledMinZero / max(ScaledMinZero);
  return(bin);
}

BinPerFile = rbindlist(lapply(unique(BinPerFile$Method), function(m) ScaleMinMax(BinPerFile[BinPerFile$Method==m, ])));

g = ggplot(BinPerFile, aes(x=BinID, y=RescaledCount, color=Method)) + geom_line();
g = g + ylab(label="Normalized Read Coverage") + xlab("Normalized Distance Along Transcript 5' -> 3' (%)");
g = g + ggtitle(paste(format(length(unique(BinFlat$Gene)), big.mark=","), " Genes", sep=""));
#g;

subBinPerFile = BinPerFile[grep("10ng_jurkat_run17_001", BinPerFile$Method), ];
g = ggplot(subBinPerFile, aes(x=BinID, y=RescaledCount, color=Method)) + geom_line();
g = g + ylab(label="Normalized Read Coverage") + xlab("Normalized Distance Along Transcript 5' -> 3' (%)");
g = g + ggtitle(paste(format(length(unique(BinFlat$Gene)), big.mark=","), " Genes", sep=""));
#g;
