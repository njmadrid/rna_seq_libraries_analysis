# create a merged splice junction table for Roger & Scott
# do it on a laptop so it won't melt their laptops
library(plyr)
library(dplyr)
library(tidyr)
library(parallel)

tabfiles <- dir(path="c:/Users/Jason/Desktop/rna_prep_comparison/tabs/", pattern=".*.out.tab$", full.names = T)

foo <- lapply(tabfiles, read.table, header=F, stringsAsFactors=F)
names(foo) <- sub("^rna_","",sub("_SJ.out.tab$","",basename(tabfiles)))
junctiontbl <- bind_rows(foo, .id = "sample")
colnames(junctiontbl) <- c("sample", "chromosome","first.intron.base","last.intron.base", "strand", "intron.motif", 
                        "annotated", "unique.mapping","multi.mapping","max.overhang")

unique.mapping <- junctiontbl %>% select(sample,chromosome, first.intron.base, last.intron.base, 
                                         strand, intron.motif, annotated, unique.mapping) %>% 
                  spread(key=sample, value=unique.mapping, fill=0)

multi.mapping <- junctiontbl %>% select(sample,chromosome, first.intron.base, last.intron.base, 
                                     strand, intron.motif, annotated, multi.mapping) %>% 
                  spread(key=sample, value=multi.mapping, fill=0)

all.mapping <- junctiontbl %>% mutate(all.mapping = unique.mapping + multi.mapping) %>%
                  select(sample,chromosome, first.intron.base, last.intron.base, 
                                      strand, intron.motif, annotated, all.mapping) %>% 
                  spread(key=sample, value=all.mapping, fill=0)

# how many novel vs annotated junctions
table(all.mapping$annotated) # 670302 vs 265411
# granted a good portion of the novel will be garbage


# get a hint what a lower threshold should be. 5 reads to count an event appears reasonable
table(junctiontbl$unique.mapping)[1:20]
table(junctiontbl$multi.mapping)[1:20]
table(junctiontbl$unique.mapping+junctiontbl$multi.mapping)[1:20]

# junctions observed
total.junctions <- apply(all.mapping[,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))
annotated.junctions <- apply(all.mapping[all.mapping$annotated==1,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))
novel.junctions <- apply(all.mapping[all.mapping$annotated==0,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))

df <- data.frame(total.junctions, annotated.junctions, novel.junctions)

View(df)

write.csv(unique.mapping,"uniquemapping.csv", row.names=F)
write.csv(multi.mapping,"multimapping.csv", row.names=F)
write.csv(all.mapping,"allmapping.csv", row.names=F)
write.csv(df, "junctionsummary.csv")
